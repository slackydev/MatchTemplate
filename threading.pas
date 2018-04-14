unit Threading;
{==============================================================================]
  Copyright © 2018, Jarl Krister Holta
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
[==============================================================================}
{$mode objfpc}{$H+}
{$modeswitch advancedrecords}
{$inline on}
interface
uses
  {$IFDEF UNIX}cthreads, cmem,{$ENDIF}
  SysUtils, Classes, core;

type
  TThreadMethod = procedure(params: PParamArray);
  TThreadId = Int32;
  
  TExecThread = class(TThread)
  protected
    FMethod: TThreadMethod;
    FParams: TParamArray;
    procedure Execute; override;
  public
    Executed: Boolean;
    constructor Create();
    procedure SetMethod(Method: TThreadMethod); inline;
    procedure SetArgument(argId:Int32; arg:Pointer); inline;
  end;
  
  TThreadArray = array of record
    Thread: TExecThread;
    Available: Boolean;
    Initialized: Boolean;
  end;


  TThreadPool = class(TObject)
    FMaxThreads: SizeInt;
    FThreads: TThreadArray;

    constructor Create(MaxThreads: SizeInt);
    destructor Free();

    function GetAvailableThread(): TThreadId;
    function  NewThread(method: TThreadMethod): TThreadId;
    procedure SetArgument(t:TThreadId; argId:Int32; arg:Pointer); inline;
    procedure SetArguments(t:TThreadId; args:array of Pointer);
    procedure Start(t:TThreadId);
    function  Executed(t:TThreadId): Boolean; inline;
    procedure Release(t:TThreadId);

    procedure MatrixFunc(Method:TThreadMethod; Args: array of Pointer; W,H:Int32; nThreads:UInt8; fallback:Int32=64*64);
    procedure TestMatrixFunc(Method:TThreadMethod; Args: array of Pointer; W,H:Int32; nThreads:UInt8);
  end;

function GetSystemThreadCount(): Int32;

//--------------------------------------------------
implementation

{$IFDEF LINUX}{$linklib c}{$ENDIF}

uses
  {$IF defined(WINDOWS)}Windows,
  {$ELSEIF defined(freebsd) or defined(darwin)}ctypes,sysctl,
  {$ELSEIF LINUX}ctypes,{$ENDIF}
  Math;

function GetSystemThreadCount(): Int32;
{$IF Defined(WINDOWS)}
var
  i: Integer;
  ProcessAffinityMask, SystemAffinityMask: DWORD_PTR;
  Mask: DWORD;
  SystemInfo: SYSTEM_INFO;
begin
  if GetProcessAffinityMask(GetCurrentProcess, ProcessAffinityMask, SystemAffinityMask)
  then begin
    Result := 0;
    for i := 0 to 31 do begin
      Mask := DWord(1) shl i;
      if (ProcessAffinityMask and Mask) <> 0 then
        Inc(Result);
    end;
  end else begin
    //can't get the affinity mask so we just report the total number of processors
    GetSystemInfo(SystemInfo);
    Result := SystemInfo.dwNumberOfProcessors;
  end;
end;
{$ELSEIF Defined(freebsd) OR Defined(darwin)}
var
  mib: array[0..1] of cint;
  len: cint;
  t: cint;
begin
  mib[0] := CTL_HW;
  mib[1] := HW_NCPU;
  len := SizeOf(t);
  fpsysctl(pchar(@mib), 2, @t, @len, Nil, 0);
  Result:=t;
end;
{$ELSEIF Defined(linux)}
function SysConf(i: cint): clong; cdecl; external name 'sysconf';
begin
  Result := SysConf(83);
end;
{$ELSE}
begin
  Result := 1;
end;
{$ENDIF}


(*----| WorkThread |----------------------------------------------------------*)
constructor TExecThread.Create();
begin
  FreeOnTerminate := True;
  Executed := False;
  inherited Create(True, 1*1024*1024);
end;

procedure TExecThread.Execute;
begin
  FMethod(@FParams);
  Executed := True;
end;

procedure TExecThread.SetMethod(Method: TThreadMethod);
begin
  FMethod := Method;
end;

procedure TExecThread.SetArgument(argId:Int32; arg:Pointer);
begin
  Self.FParams[argId] := arg;
end;


(*----| ThreadPool |----------------------------------------------------------*)
constructor TThreadPool.Create(MaxThreads:SizeInt);
var i:Int32;
begin
  Self.FMaxThreads := MaxThreads;
  SetLength(self.FThreads, FMaxThreads);
  for i:=0 to High(self.FThreads) do
  begin
    self.FThreads[i].Available   := True;
    self.FThreads[i].Initialized := False;
  end;
end;

destructor TThreadPool.Free();
var i:Int32;
begin
  for i:=0 to High(self.FThreads) do
    self.Release(i);
end;

function TThreadPool.GetAvailableThread(): TThreadId;
var i:Int32;
begin
  for i:=0 to High(self.FThreads) do
    if self.FThreads[i].Available then
      Exit(TThreadId(i));
  raise Exception.Create('TThreadPool.GetAvailableThread: No free execution threads'); 
end;

function TThreadPool.NewThread(method: TThreadMethod): TThreadId;
begin
  result := GetAvailableThread();
  if self.FThreads[result].Thread <> nil then
    Self.Release(result);

  self.FThreads[result].Thread := TExecThread.Create();
  self.FThreads[result].Available := False;
  self.FThreads[result].Thread.SetMethod(method);
end;

procedure TThreadPool.SetArgument(t:TThreadId; argid:Int32; arg:Pointer);
begin
  self.FThreads[t].Thread.SetArgument(argid, arg);
end;

procedure TThreadPool.SetArguments(t:TThreadId; args:array of Pointer);
var arg:Int32;
begin
  for arg:=0 to High(args) do
    self.FThreads[t].Thread.SetArgument(arg, args[arg]);
end;

procedure TThreadPool.Start(t:TThreadId);
begin
  self.FThreads[t].Thread.Start;
end;

function TThreadPool.Executed(t:TThreadId): Boolean;
begin
  if (self.FThreads[t].Thread = nil) or (self.FThreads[t].Available) then
     Result := False
  else
    Result := self.FThreads[t].Thread.Executed = True;
end;

procedure TThreadPool.Release(t:TThreadId);
begin
  self.FThreads[t].Available := True;
  if self.FThreads[t].Thread <> nil then
  begin
    self.FThreads[t].Thread.Terminate();
    self.FThreads[t].Thread := nil;
  end;
end;

(*----| Functions |----------------------------------------------------------*)
procedure TThreadPool.MatrixFunc(Method:TThreadMethod; Args: array of Pointer; W,H:Int32; nThreads:UInt8; fallback:Int32=64*64);
var
  i,lo,hi,step: Int32;
  thread: array of record id: TThreadId; box: TBox; end;
  params: TParamArray;
  area: TBox;
begin
  if (W*H < fallback) or (nThreads=1) then
  begin
    area := Box(0,0,W-1,H-1);
    params := args;
    params[length(args)] := @area;
    Method(@params);
    Exit();
  end;

  nThreads := Max(1, nThreads);
  SetLength(thread, nThreads);
  
  lo := 0;
  step := (H-1) div nThreads; 
  for i:=0 to nThreads-1 do
  begin
    hi := Min(H-1, lo + step);
    
    thread[i].box := Box(0, lo, w-1, hi);
    thread[i].id := Self.NewThread(Method);
    Self.SetArguments(thread[i].id, Args);
    Self.SetArgument(thread[i].id, length(args), @thread[i].box);
    Self.Start(thread[i].id);
	
    if hi = H-1 then
    begin
      nThreads := i+1;
      Break;
    end;
    lo := hi + 1;
  end;

  for i:=0 to nThreads-1 do
  begin
    while not Self.Executed(thread[i].id) do Sleep(0);
    Self.Release(thread[i].id);
  end;
end;

procedure TThreadPool.TestMatrixFunc(Method:TThreadMethod; Args: array of Pointer; W,H:Int32; nThreads:UInt8);
var
  i,lo,hi,step: Int32;
  params: TParamArray;
  area: TBox;
begin
  nThreads := Max(1, nThreads);
  lo := 0;
  step := (H-1) div nThreads;
  for i:=0 to nThreads-1 do
  begin
    hi := Min(H-1, lo + step);

    area := Box(0, lo, w-1, hi);
    params := args;
    params[length(args)] := @area;
    Method(@params);

    if hi = H-1 then Break;
    lo := hi + 1;
  end;
end;


end.
