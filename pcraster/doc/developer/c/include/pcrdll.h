#ifndef  INCLUDED_PCRDLL
#define  INCLUDED_PCRDLL

/*!
   \file
    defines to tag windows dll exported funcs

   A Dll should export all functions callable from any language.
   Visual Basic only supports stdcall. Java does stdcall on Windows and cdecl
   on UNIX.

   In addition to the macro PCR_DLL_FUNC a definition should be made
   in the .def file when linking for windows, to make the function
   visible as is, without decoration (e.g. prefix _ and @4 suffix).
   This not done yet, which is only a problem for VB6 that needs to use the
   decorated functions as is:
   Private Declare Function ScriptError Lib "pcrme.dll" Alias "_pcr_ScriptError@4" (ByVal lngScriptHandle As Long) As Long
   (For an example see libs/pcrme/pcrme.def)
 */

#if defined(WIN32) || defined(_WIN32)
/* seems to work like this in MSVC/MINGW/BORLANDC
 */
/* # define PCR_DLL_FUNC(retType) __declspec(dllexport) retType __stdcall */
# define PCR_DLL_FUNC(...) __declspec(dllexport) __VA_ARGS__ __stdcall
# define PCR_DLL_C             __declspec(dllexport)
# define PCR_DLL_CLASS         __declspec(dllexport)
#else
/* # define PCR_DLL_FUNC(retType) retType */
# define PCR_DLL_FUNC(...) __VA_ARGS__
# define PCR_DLL_C
# define PCR_DLL_CLASS
#endif

#endif
