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

  --

  Standalone unit for basic parallel processing of arrays (needs cpuinfo).
[==============================================================================}
{$mode objfpc}{$H+}
{$modeswitch advancedrecords}
{$inline on}
interface
uses
  {$ifdef unix}cthreads,cmem,{$endif} // c memory manager can be notably faster on some platforms
  SysUtils, Classes;

type
  PParamArray = ^TParamArray;
  TParamArray = array[Word] of Pointer;

  TThreadMethod = procedure(Params: PParamArray; iLow, iHigh: Int32);
  TThreadId = Int32;

  EThreadState = (tsSleeping, tsAwaiting, tsWorking, tsReady);

  TExecThread = class(TThread)
  private
    FMethod: TThreadMethod;
    FParams: TParamArray;
    FRangeLo, FRangeHi: Int32;
    procedure Execute; override; // hide this method to avoid confusion
  public
    State: EThreadState;

    constructor Create();
    procedure SetMethod(Method: TThreadMethod);
    procedure SetArgument(argId:Int32; arg:Pointer);
    procedure SetArguments(args: array of Pointer);
    function Run(iLow, iHigh: Int32): Boolean;
    function Completed(): Boolean; inline;
  end;
  TThreadArray = array of TExecThread;

  TThreadPool = class(TObject)
  private
    FMaxThreads: SizeInt;
    FThreads: TThreadArray;
    function  FindAvailableThread(): TThreadId; inline;
    function  PrepareThread(Method: TThreadMethod): TExecThread; inline;
  public  
    constructor Create(MaxThreads: SizeInt);
    destructor Free();
    procedure DoParallel(Method: TThreadMethod; Args: array of Pointer; iLow, iHigh: Int32; nThreads: UInt8; Fallback: Boolean=False);
  end;

 
var
  ThreadPool: TThreadPool;



implementation

uses
  Math, cpuinfo;


// ----------------------------------------------------------------------------
// Execution thread

constructor TExecThread.Create();
begin
  FreeOnTerminate := True;
  State := tsSleeping;
  inherited Create(False, 1024*512); //default = 4MiB, we set 512KiB
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

procedure TExecThread.SetArgument(argId: Int32; arg: Pointer);
begin
  Self.FParams[argId] := arg;
end;

procedure TExecThread.SetArguments(args: array of Pointer);
var arg:Int32;
begin
  for arg:=0 to High(args) do
    self.FParams[arg] := args[arg];
end;

function TExecThread.Run(iLow, iHigh: Int32): Boolean;
begin
  Result := Self.FMethod <> nil;
  if Result then
  begin
    FRangeLo := iLow;
    FRangeHi := iHigh;

    self.State := tsReady;
  end;
end;

function TExecThread.Completed(): Boolean;
begin
  Result := self.State = tsSleeping;
end;



// ----------------------------------------------------------------------------
// ThreadPool

constructor TThreadPool.Create(MaxThreads: SizeInt);
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

function TThreadPool.PrepareThread(Method: TThreadMethod): TExecThread;
begin
  Result := self.FThreads[FindAvailableThread()];
  Result.State := tsAwaiting;
  Result.SetMethod(Method);
end;

procedure TThreadPool.DoParallel(Method: TThreadMethod; Args: array of Pointer; iLow, iHigh: Int32; nThreads: UInt8; Fallback: Boolean=False);
var
  i,step,A,B: Int32;
  Thread: array of TExecThread;
  Params: TParamArray;
begin
  if (fallback) or (nThreads=1) then
  begin
    Params := Args;
    Method(@Params, iLow, iHigh);
    Exit();
  end;

  nThreads := Max(1, nThreads);
  SetLength(thread, nThreads);

  A := iLow;
  step := Max(1, (iHigh+1) div nThreads);

  for i:=0 to nThreads-1 do
  begin
    B := Min(iHigh, A + step);

    Thread[i] := Self.PrepareThread(Method);
    Thread[i].SetArguments(Args);
    Thread[i].Run(A,B);

    if B = iHigh then
    begin
      nThreads := i+1;
      Break;
    end;
    A := B + 1;
  end;

  for i:=0 to nThreads-1 do
    while not Thread[i].Completed do
      Sleep(0);
end;


initialization
  WriteLn('>>> Spawning ',GetSystemThreadCount(),' worker threads');
  ThreadPool := TThreadPool.Create(GetSystemThreadCount());

end.
