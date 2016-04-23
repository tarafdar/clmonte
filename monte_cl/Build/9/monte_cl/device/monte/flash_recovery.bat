@rem ACL utility command
@rem 
@rem Example:  acl help

@echo off
@if not exist "%WINDIR%\system32\chcp.com" goto after_chcp
for /f "tokens=2 delims=:." %%x in ('chcp') do set saved_chcp=%%x
@rem Change console code to US english
@chcp 437>nul
:after_chcp
@echo on
@if not "%ALTERAOCLSDKROOT%" == "" goto check_bin
@echo Error: You must set environment variable ALTERAOCLSDKROOT to the absolute path of root of the Altera SDK for OpenCL software installation
@set ERRORLEVEL=1
@goto end
:check_bin
@if not "%QUARTUS_ROOTDIR%" == "" goto check_qtcl
@echo Error: You must set environment variable QUARTUS_ROOTDIR to the absolute path of the quartus subdirectory inside ACDS
@set ERRORLEVEL=2
@goto end
:check_qtcl
@if exist "%QUARTUS_ROOTDIR%/common/tcl/internal" goto preamble_ok
@echo Error: You must set environment variable QUARTUS_ROOTDIR to the absolute path of the quartus subdirectory inside ACDS
@set ERRORLEVEL=2
@goto end
:preamble_ok
@set PERL5LIB=%ALTERAOCLSDKROOT%\share\lib\perl;%ALTERAOCLSDKROOT%\share\lib\perl\5.8.8;%PERL5LIB%
@echo About to call perl.
%QUARTUS_ROOTDIR%\bin64\perl\bin\perl.exe flash_recovery.pl base.aocx
@echo after call to perl
@if not exist "%WINDIR%\system32\chcp.com" goto end
@rem Restore console settings
@chcp %saved_chcp%>nul
:end
