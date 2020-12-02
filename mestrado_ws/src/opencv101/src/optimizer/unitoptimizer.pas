unit unitoptimizer;

{$link optimizer.o}
{$link meschach.a}
{$linklib c}
{$linklib m}
{$L donlp2.o}
{$L newx.o}
{$L user_eval.o}
{$mode objfpc}{$H+}


interface

uses CTypes;


procedure optimizer(Params:Array Of Ctypes.cdouble; x00: Array Of Ctypes.cdouble; u00: Array Of Ctypes.cdouble; TC1: Array Of Ctypes.cdouble; TC2: Array Of Ctypes.cdouble; var result: Array of Ctypes.cdouble); cdecl; external;


implementation

end.
