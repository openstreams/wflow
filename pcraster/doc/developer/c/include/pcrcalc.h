#ifndef INCLUDED_PCR_PCRCALC
#define INCLUDED_PCR_PCRCALC

#ifndef INCLUDED_PCRDLL
#include "pcrdll.h"
#define INCLUDED_PCRDLL
#endif

/*
 * All rights reserved
 *  by PCRaster Environmental Software
 */

#ifndef __cplusplus
// * then generate only the external docs
//    by make apiDocs
// * Docs are in header not in body: header is distributed in SDK,
//   body not
// * NOTE doxygen wants PCR_DLL_FUNC macro on single line
// * Check calc_clientinterfacetest.cc foo() for
//    verbatim piece to copy into user docs

/*! \brief handle to a script returned by pcr_createScriptFromTextFile(), pcr_createScriptFromTextString(), pcr_createScriptFromXMLFile() or pcr_createScriptFromXMLString()
 *
 * PcrScript* is an opaque pointer to the runtime structure
 * associated with a single script.
 * The structure should be deleted with pcr_destroyScript()
 */
 typedef struct PcrScript { char dummy; } PcrScript;
#else
// attach to Pcrcalc class in PCRasterModelEngine.dll
 typedef struct Pcrcalc PcrScript;
 namespace calc {
   class ASTScript;
 }
 extern "C" {
#endif

/*!
   Both stepwise functions, pcr_ScriptExecuteInitialStepMemory() and pcr_ScriptExecuteNextTimeStepMemory(),
   have a dataTransferArray argument used for transferring data.
   dataTransferArray is a low level interface. The memory layout of DataTransferArray is specified
   in the script XML document.

   dataTransferArray is the address of an array of addresses (in C terms:
   an array of void pointers). The index in the array is the nominal value
   presented  in the definition/{scriptInput/scriptOutput}/memoryExchange element.
   The array must have an entry for each memoryExchange element present.
   If the entry is not 0, then the entry is the address of a transfer buffer.

   The semantics of the transfer buffer depends on if the memoryExchange
   element appears in a scriptInput or scriptOutput element:
    - <b>scriptInput</b>: The buffer does contain a value for the item, when the stepwise function is called. The buffer is not modified when the stepwise function returns. An entry may be 0 if the corresponding memoryExchange element is not used in either the initial or dynamic section for which the stepwise function is called. If it used, but 0 is passed a runtime error will occur (see pcr_ScriptError()).
    - <b>scriptOutput</b>: If the buffer is not 0 when the stepwise function is called then the caller has reserved enough space at the buffer address. The stepwise function will copy data into the buffer. If the buffer is 0 the stepwise function will allocate it (see <i>allocation by the function</i> below).
    .


    \section functionAllocation Allocation by the Function
   If the buffer is 0 the stepwise function will allocate it. The buffer allocated should not be deleted/freed by the caller nor should the address by passed back in a successive call to a stepwise function. A successive call to a stepwise function or pcr_destroyScript() will delete the buffer. Thus the buffer is only valid until such a next call.

   <b>Rule</b>: Buffers allocated by the caller should be deleted by the caller. Buffers allocated by the stepwise functions will be deleted by a next stepwise function or pcr_destroyScript() call.

    \section transferBuffer Transfer Buffer layout

    The transfer buffer is either:
     - a field buffer, if the definition element the memoryExchange element is appearing in holds a field element.
     - an indexedArray buffer, if the definition element the memoryExchange element is appearing in holds an indexedArray element.
     - a textString buffer if the definition element's name is used in a textStatistics element.
     - a timeoutputRow buffer if the definition element's name is assigned the result of a timeoutput function in the model.

    \subsection fieldBufferMemoryLayout Field Buffer Memory Layout

    A field buffer consist of one the following elementary types:
    - unsigned char (UINT1): 8-bit byte for Boolean and Ldd.
    - int (INT4): 32 bit integer for Nominal and Ordinal.
    - float (REAL4): 32 bit float for Scalar and Directional.
    .

    The field buffer is either:
     - a single value if the definition element contains the spatialType element with the value NonSpatial
     - an array of size areaMap/nrRows * areaMap/nrCols if the
       field definition element contains the spatialType element with the
       value Spatial
     .

    \subsection indexedArrayBufferMemoryLayout IndexedArray Buffer Memory Layout

    An indexedArray buffer consists of two parts: a header followed by the data.

    The header consists of 1+NDIM values of the type specified in the indexedArray/dimensionDataType element. The first value of the header represent the number of dimensions (NDIM). Successive values will give the length of each dimension.

    The data consists of X values of the type specified in the indexedArray/valueDataType. X is the multiplication of all dimension lengths.

    Currently only a 1 dimensional array is needed to use the lookup... functions with an indexedArray element. For example an array of 1 by 3 with the array values 1.33, 2.33 and 3.33 will have the following layout:
    \verbatim
      1    # number of dimensions
      3    # length of dimension 1
      1.33
      2.33
      3.33
    \endverbatim

    \subsection textStringMemoryLayout textString Buffer Memory Layout

    A 0-terminated ANSI character string. Since the buffer is of unknown length upon calling it is recommended to let the function allocate this buffer (pass 0 in, get address out).

    \subsection timeoutputRowLayout timeoutput Buffer Memory Layout

    In essence this a 1 dimensional array of values as they would appear on a row of a pcrcalc timeoutput tss file, with the time column removed.

    The memory buffer does however prepend a header of 4 unsigned integers (INT4):
    - 0. Always holds value 1 (to allow other versions in the future).
    - 1. Constant (matching the CSF_CR enum of the PCRaster file format API)
       identifying the type of the row entries:
       -  0 (0x00) for UINT1 (unsigned char)
       - 38 (0x26) for INT4 (signed int)
       - 90 (0x5A) for REAL4 (float)
       .
    - 2. Always holds value 1 (the number of dimensions of this buffer).
    - 3. The number of entries following (length of dimension 1).

    After the header the row entries will follow which are of the type identified by header item 1 and the number of entries by header item 3.

    Since the buffer is of unknown length upon calling it is recommended to let the function allocate this buffer (pass 0 in, get address out).

 */
typedef void *DataTransferArray[];


//! create the PcrScript object by supplying a filename of a text script file
/*!
 *  pcr_createScriptFromTextFile allocates some internal memory to set up error handling.
 *  No other actions, such as checking if the fileName exists or the script has a valid syntax,
 *  are taken here. All these actions are deferred
 *  to a next call that needs the script such as pcr_ScriptExecute() 
 *  or pcr_ScriptXMLReflection().
 *
 * \param fileName  C-string (Ansi 1-byte char set) with name of script file.
 *
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(PcrScript*)
#else
 PcrScript*
#endif
  pcr_createScriptFromTextFile(const char* fileName);

//! create the PcrScript object by supplying a filename of a XML file with the script root element
/*!
 *  pcr_createScriptFromXMLFile allocates some internal memory to set up error handling.
 *  No other actions, such as checking if the file exists or the script has a valid syntax,
 *  are taken here. All these actions are deferred
 *  to a next call that needs the script such as pcr_ScriptExecute() 
 *  or pcr_ScriptXMLReflection().
 *
 * \param fileName  C-string (Ansi 1-byte char set) with name of script file
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(PcrScript*)
#else
 PcrScript*
#endif
  pcr_createScriptFromXMLFile(const char* fileName);

//! create the PcrScript object by supplying the textual script as string
/*!
 *  pcr_createScriptFromTextFile allocates some internal memory to set up error handling.
 *  No other actions, such as checking if the script has a valid syntax,
 *  are taken here. Al these actions are deferred
 *  to a next call that needs the script such as pcr_ScriptExecute() 
 *  or pcr_ScriptXMLReflection().
 *
 * \param scriptContents  C-string (Ansi 1-byte char set) with script contents
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(PcrScript*)
#else
 PcrScript*
#endif
  pcr_createScriptFromTextString(const char* scriptContents);

//! create the PcrScript object by supplying the XML as a string with the script root element
/*!
 *  pcr_createScriptFromXMLString allocates some internal memory to set up error handling.
 *  No other actions, such as checking if the script has a valid syntax,
 *  are taken here. Al these actions are deferred
 *  to a next call that needs the script such as pcr_ScriptExecute() 
 *  or pcr_ScriptXMLReflection().
 *
 * \param xmlString  C-string (Ansi 1-byte char set) with script contents
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(PcrScript*)
#else
 PcrScript*
#endif
  pcr_createScriptFromXMLString(const char* xmlString);

//! clean up PcrScript object
/*!
 *  cleans up all dynamic data kept for \a script
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(void)
#else
 void
#endif
  pcr_destroyScript(PcrScript *script);

/*! check if an error has occured
 * \param script  the script object
 *
 * Read \ref errorHandling on how errors are handled.
 *
 * \returns  0 if no error present for \a script, otherwise the
 *             size of the string (excl. terminating 0) returned by
 *             pcr_ScriptErrorMessage()
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(int)
#else
 int
#endif
  pcr_ScriptError(PcrScript *script);

/*! get error message
 *  \param script  the script object
 *
 * Read \ref errorHandling on how errors are handled.
 *
 * \returns a C-string (Ansi 1-byte char set), empty string if no error present for \a script, otherwise the error message.
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(const char*)
#else
 const char*
#endif
  pcr_ScriptErrorMessage(PcrScript *script);


//! PCRaster script XML document
/*!
 * \param script         the script object
 *
 * \returns 0 in case of error, C-string (Ansi 1-byte char set) with document otherwise.
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(const char*)
#else
 const char*
#endif
  pcr_ScriptXMLReflection(PcrScript *script);

//!  Execute entire script
/*!
 * This execute variant expects all data to present as files on the filesystem, in contrast
 * to the pcr_ScriptExecuteInitialStepMemory() and pcr_ScriptExecuteNextTimeStepMemory() setup.
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(void)
#else
 void
#endif
  pcr_ScriptExecute(PcrScript *script);

//! Execute the initial section only
/*!
  \param script  the script object
  \param dataTransferArray see  DataTransferArray

   If the model is a static model, the entire model is executed. Repetive calls will
   do nothing; for each script object the initial section can only be called once.

   Calling this function more than once will result in an error.

   \returns
    - -1  model terminated in error, check pcr_ScriptErrorMessage(PcrScript *script)
    -  0  initial step is executed, model finished thus it is a static model.
    -  1  initial step is executed, model not yet finished, thus it is a dynamic model.
    .

  \warning
    The function will not put an output into \a dataTransferArray if that output is also
    re-computed in the dynamic section:
    \verbatim
    initial
     # pcr_ScriptExecuteInitialStepMemory will *NOT* put value of out in dataTransferArray
     out = in * 1;
    dynamic
     # pcr_ScriptExecuteNextTimeStepMemory will put value of out in dataTransferArray
     out =  f(out, ...);
    \endverbatim
    This is consistent how pcrcalc deals with output map stacks versus single map output.
    If you need the initial value use an auxiliary variable (no performance penalty):
    \verbatim
    initial
     # pcr_ScriptExecuteInitialStepMemory will *NOT* put value of out in dataTransferArray
     out = in * 1;
     # pcr_ScriptExecuteInitialStepMemory cat put this value in dataTransferArray
     staticOut = out;
    dynamic
     # pcr_ScriptExecuteNextTimeStepMemory will put value of out in dataTransferArray
     out =  f(out, ...);
    \endverbatim
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(int)
#else
  int
#endif
  pcr_ScriptExecuteInitialStepMemory(
    PcrScript *script,
    DataTransferArray dataTransferArray);

//! Execute a timestep
/*!
   \param script  the script object
   \param dataTransferArray

   Execute next time step.
   Calling this function without first calling pcr_ScriptExecuteInitialStepMemory() will result in an error.

   \returns
     - -1  model terminated in error, check pcr_ScriptErrorMessage(PcrScript *script)
     -  0  a timestep is executed, model finished; this call executed the last timestep.
     -  1  a timestep is executed, model not yet finished, more steps may follow.
     .
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(int)
#else
 int
#endif
  pcr_ScriptExecuteNextTimeStepMemory(PcrScript *script,
                                      DataTransferArray dataTransferArray);


//! End execution started with pcr_ScriptExecuteInitialStepMemory() and pcr_ScriptExecuteNextTimeStepMemory() sequence.
/*!
 * \param script  the script object
 *
 * \returns
 *   - -1  finishing raised an error, check pcr_ScriptErrorMessage(PcrScript *script)
 *   -  0  finished with no error
 *   .
 *
 * pcr_destroyScript() will perform the same actions but pcr_destroyScript() can not return an error message.
 * Only needed if caller wants to terminate before end of model or if model does
 * not know the number of timesteps itself.
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(int)
#else
 int
#endif
  pcr_ScriptExecuteFinish(PcrScript *script);

/* NOT EXPOSED
  \brief Release all memory that has been allocated by any of the execute functions that involve a DataTransferArray with output elements requesting memory allocation by the execute function.

   \param script  the script object

   This function does nothing if no memory allocations have been requested.
   Memory allocations that have not been release are also released when
   pcr_destroyScript() is called.

   \returns
     - -1  bad \a script pointer
     -  0  no error
     .

#ifndef DOXYGEN_SHOULD_SKIP_THIS
 PCR_DLL_FUNC(int)
#else
 int
#endif
  pcr_ScriptReleaseAllAllocatedMemory(PcrScript *script);
*/

#ifdef __cplusplus
  //! not part of the API
  calc::ASTScript const& pcr_internalScript(PcrScript *script);
 }
#endif

#endif
