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
  SysUtils, Classes, core;

type
  TThreadMethod = procedure(Params: PParamArray; iLow, iHigh: Int32);
  TThreadId = Int32;

  EThreadState = (tsSleeping, tsAwaiting, tsWorking, tsReady);

  TExecThread = class(TThread)
  protected
    FMethod: TThreadMethod;
    FParams: TParamArray;
    FRangeLo, FRangeHi: Int32;

    procedure Execute; override;
  public
    State: EThreadState;

    constructor Create();
    procedure SetMethod(Method: TThreadMethod); inline;
    procedure SetArgument(argId:Int32; arg:Pointer); inline;
    procedure SetRange(Low, High: Int32); inline;
  end;
  
  TThreadArray = array of TExecThread;

  TThreadPool = class(TObject)
    FMaxThreads: SizeInt;
    FThreads: TThreadArray;

    constructor Create(MaxThreads: SizeInt);
    destructor Free();

    function  FindAvailableThread(): TThreadId; inline;
    function  NewThread(method: TThreadMethod): TThreadId; inline;
    procedure SetArgument(t: TThreadId; argId: Int32; arg: Pointer); inline;
    procedure SetArguments(t: TThreadId; args: array of Pointer); {inline;}
    procedure SetRange(t: TThreadId; Low, High: Int32); inline;

    procedure Start(t: TThreadId);
    function  Executed(t: TThreadId): Boolean; inline;

    procedure DoParallel(Method: TThreadMethod; Args: array of Pointer; iLow, iHigh: Int32; nThreads: UInt8; fallback: Boolean=False);
  end;

var
  ThreadPool: TThreadPool;

//--------------------------------------------------
implementation

uses
  Math;


(*----| WorkThread |----------------------------------------------------------*)
constructor TExecThread.Create();
begin
  FreeOnTerminate := True;
  State := tsSleeping;
  inherited Create(False, 1024*1024);
end;

procedure TExecThread.Execute();
label startloc;
begin
  startloc:

  Self.Priority := tpIdle;
  while (Self.State <> tsReady) do Sleep(1);
  Self.Priority := tpNormal;

  State := tsWorking;
  FMethod(@FParams, FRangeLo, FRangeHi);
  State := tsSleeping;

  goto startloc;
end;

procedure TExecThread.SetMethod(Method: TThreadMethod);
begin
  FMethod := Method;
end;

procedure TExecThread.SetArgument(argId:Int32; arg:Pointer);
begin
  Self.FParams[argId] := arg;
end;

procedure TExecThread.SetRange(Low, High: Int32);
begin
  FRangeLo := Low;
  FRangeHi := High;
end;


(*----| ThreadPool |----------------------------------------------------------*)
constructor TThreadPool.Create(MaxThreads:SizeInt);
var i:Int32;
begin
  Self.FMaxThreads := MaxThreads;
  SetLength(self.FThreads, FMaxThreads);
  for i:=0 to High(self.FThreads) do
    self.FThreads[i] := TExecThread.Create();
end;

destructor TThreadPool.Free();
var i:Int32;
begin
  for i:=0 to High(self.FThreads) do
  begin
    if self.FThreads[i] <> nil then
    begin
      self.FThreads[i].Terminate();
      self.FThreads[i] := nil;
    end;
  end;
end;

function TThreadPool.FindAvailableThread(): TThreadId;
var i:Int32;
begin
  while True do
  begin
    for i:=0 to High(self.FThreads) do
      if (self.FThreads[i].State = tsSleeping) then
        Exit(TThreadId(i));
    Sleep(0);
 end;
end;

function TThreadPool.NewThread(method: TThreadMethod): TThreadId;
begin
  Result := FindAvailableThread();
  self.FThreads[result].State := tsAwaiting;
  self.FThreads[result].SetMethod(method);
end;

procedure TThreadPool.SetArgument(t: TThreadId; argid: Int32; arg: Pointer);
begin
  self.FThreads[t].SetArgument(argid, arg);
end;

procedure TThreadPool.SetArguments(t: TThreadId; args: array of Pointer);
var arg:Int32;
begin
  for arg:=0 to High(args) do
    self.FThreads[t].SetArgument(arg, args[arg]);
end;

procedure TThreadPool.SetRange(t: TThreadId; Low, High: Int32);
begin
  self.FThreads[t].SetRange(Low, High);
end;

procedure TThreadPool.Start(t: TThreadId);
begin
  self.FThreads[t].State := tsReady;
end;

function TThreadPool.Executed(t: TThreadId): Boolean;
begin
  Result := self.FThreads[t].State = tsSleeping;
end;


// -----------------------------------------------------------------------------
// Functions

procedure TThreadPool.DoParallel(Method:TThreadMethod; Args: array of Pointer; iLow, iHigh: Int32; nThreads: UInt8; fallback: Boolean=False);
var
  i,step,A,B: Int32;
  thread: array of TThreadId;
  params: TParamArray;
begin
  if (fallback) or (nThreads=1) then
  begin
    params := Args;
    Method(@params, iLow, iHigh);
    Exit();
  end;

  nThreads := Max(1, nThreads);
  SetLength(thread, nThreads);
  
  A := iLow;
  step := Max(1, (iHigh+1) div nThreads);

  for i:=0 to nThreads-1 do
  begin
    B := Min(iHigh, A + step);

    thread[i] := Self.NewThread(Method);
    Self.SetArguments(thread[i], Args);
    Self.SetRange(thread[i], A,B);
    Self.Start(thread[i]);

    if B = iHigh then
    begin
      nThreads := i+1;
      Break;
    end;
    A := B + 1;
  end;

  for i:=0 to nThreads-1 do
    while not Self.Executed(thread[i]) do
      Sleep(0);
end;


initialization
  ThreadPool := TThreadPool.Create(64);

end.
