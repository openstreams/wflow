using System;
using NUnit.Framework;


namespace pcr.Test {


/// <summary>
///  Test class to demonstrates the C-API
/// </summary>
[TestFixture]
public class APITest : API {

  // 0-ptr data, not ready to run
  IntPtr[] data0      = { IntPtr.Zero, IntPtr.Zero };

  [Test] public void TestCreate()
  {
   IntPtr s=pcr_createScriptFromXMLFile("memoryOnlyIO_7.xml");
   Assert.IsNotNull(s);
   pcr_destroyScript(s);
  }

  [Test] public void TestRuntimeError()
  { // Runtime error passing 0-ptrs
   IntPtr s=pcr_createScriptFromXMLFile("memoryOnlyIO_7.xml");
   Assert.IsTrue(pcr_ScriptError(s)==0);
   pcr_ScriptExecuteInitialStepMemory(s, data0);

   Assert.IsFalse(pcr_ScriptError(s)==0);
   Assert.IsTrue(-1 != pcr_ScriptErrorMessage(s).IndexOf(
              "memInput: 0-ptr data buffer passed"));

   pcr_destroyScript(s);
  }

  [Test] unsafe public void TestExecution_memoryOnlyIO_7()
  { // Thunderbirds are go!

   IntPtr s=pcr_createScriptFromXMLFile("memoryOnlyIO_7.xml");
   Assert.IsTrue(pcr_ScriptError(s)==0);


   float *input  = stackalloc float[25];
   float *output = stackalloc float[25];
   for(int i=0; i < 25; i++) {
     input[i] =4.5F;
     // no need for init of output
     // except to make sure test is working correct
     output[i]=0F;
   }
   IntPtr[] data = { (IntPtr)null, (IntPtr)null};
   data[0]  = (IntPtr)input;
   data[1] = (IntPtr)output;


   pcr_ScriptExecuteInitialStepMemory(s, data);
   Assert.IsTrue(pcr_ScriptError(s)==0);

   // check if all are 4.5 + 3
   for(int i=0; i < 25; i++)
     Assert.AreEqual(7.5F,output[i]);

   Assert.IsTrue(pcr_ScriptError(s)==0);

   Assert.IsNotNull(s);
   pcr_destroyScript(s);
  }
}

} // namespace pcr.Test
