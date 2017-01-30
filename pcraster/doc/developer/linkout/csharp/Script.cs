using System.Runtime.InteropServices; // DllImport
using System; // IntPtr, String

using System.IO;
using System.Diagnostics;
using System.Xml;
using System.Xml.Serialization;


namespace pcr {

public class API {


    // The PCR C Application Programming Interface

    [DllImport("PCRasterModelEngine.dll",
    CharSet=CharSet.Ansi)]
    public static extern IntPtr pcr_createScriptFromTextFile(String scriptName);

    [DllImport("PCRasterModelEngine.dll",
    CharSet=CharSet.Ansi)]
    public static extern IntPtr pcr_createScriptFromTextString(String script);

    [DllImport("PCRasterModelEngine.dll",
    CharSet=CharSet.Ansi)]
    public static extern IntPtr pcr_createScriptFromXMLFile(String scriptName);

    [DllImport("PCRasterModelEngine.dll",
    CharSet=CharSet.Ansi)]
    public static extern IntPtr pcr_createScriptFromXMLString(String script);

    [DllImport("PCRasterModelEngine.dll",
    CharSet=CharSet.Ansi)]
    public static extern void pcr_destroyScript(IntPtr script);


    public static String pcr_ScriptXMLReflection(IntPtr script)
     {
       // rewrap for string conversion
       IntPtr s = pcr_ScriptXMLReflectionDll(script);
       String xmlString=Marshal.PtrToStringAnsi(s);
       // perform optional validation
       verifyXml(xmlString);
       return xmlString;
     }

    [DllImport("PCRasterModelEngine.dll",
     CharSet=CharSet.Ansi)]
     public static extern Int32  pcr_ScriptError(IntPtr script);

     public static String pcr_ScriptErrorMessage(IntPtr script) {
       // rewrap for string conversion
       IntPtr s = pcr_ScriptErrorMessageDll(script);
       return Marshal.PtrToStringAnsi(s);
     }



    [DllImport("PCRasterModelEngine.dll",
    CharSet=CharSet.Ansi)]
    public static extern void pcr_ScriptExecuteInitialStepMemory(
       IntPtr script, IntPtr[] data);

    [DllImport("PCRasterModelEngine.dll",
     CharSet=CharSet.Ansi,
     EntryPoint="pcr_ScriptXMLReflection")]
     private static extern IntPtr pcr_ScriptXMLReflectionDll(IntPtr script);

    [DllImport("PCRasterModelEngine.dll",
     CharSet=CharSet.Ansi,
     EntryPoint="pcr_ScriptErrorMessage")]
     private static extern IntPtr pcr_ScriptErrorMessageDll(IntPtr script);




    // some support/internal logic for API

    // create pcrxml::Script according to XML-schema from string
    public static pcrxml.Script xmlStringToXMLScript(String xmlString) {
      // Create an instance of the XmlSerializer specifying type and namespace.
      XmlSerializer s  = new XmlSerializer(typeof(pcrxml.Script));
      StringReader  sr = new StringReader(xmlString);
      XmlTextReader r  = new XmlTextReader(sr);
      // Declare an object variable of the type to be deserialized.
      // Use the Deserialize method to restore the object's state.
      return (pcrxml.Script)s.Deserialize(r);
    }
    // create pcrxml::Script according to XML-schema from file
    public static pcrxml.Script xmlFileToXMLScript(String xmlFileName) {
      // Create an instance of the XmlSerializer specifying type and namespace.
      XmlSerializer s  = new XmlSerializer(typeof(pcrxml.Script));
      StreamReader  sr = new StreamReader(xmlFileName);
      XmlTextReader r  = new XmlTextReader(sr);
      // Declare an object variable of the type to be deserialized.
      // Use the Deserialize method to restore the object's state.
      return (pcrxml.Script)s.Deserialize(r);
    }

    private static bool d_doXmlValidation=true;

    public static bool doXmlValidation {
      get { return d_doXmlValidation; }
      set { d_doXmlValidation=value; }
    }

    public  static bool validXml(String xmlString) {
      verifyXml(xmlString);
      return true;
    }
    [Conditional("DEBUG")]
    public  static void verifyXml(String xmlString) {
      if (doXmlValidation) {
        XsdValidation xv=new XsdValidation();
        Debug.Assert(xv.validateString(xmlString));
      }
    }

}


/// <summary>
///  support class for pcr_ScriptXMLReflection API
/// </summary>
public class XMLReflection {
   private pcrxml.Script d_script;

   public XMLReflection(String xmlString) {
     d_script=API.xmlStringToXMLScript(xmlString);
   }

   public pcrxml.Script script {
     get {

       return d_script;
     }
   }
}

/*
 * 
 * /// <summary>
 * /// Tiny engine capable of running model scripts. Convienence class that will wrap the non object oriented C API into a class. 
 * //// Current functionality is minimal, in development for OpenMI and ArcGIS 8+ interface.
 * /// </summary>
 * public class Script : API {
 * 
 *     public class Exception : System.Exception {
 *       public Exception(IntPtr script):
 *         base("Script.Exception:\n"+pcr_ScriptErrorMessage(script)) {
 *       }
 *     }
 * 
 *     private IntPtr  d_script = IntPtr.Zero;
 * 
 *     public Script(String scriptFileOrString, bool file) {
 *       if (file)
 *         d_script=pcr_createScriptFromFile(scriptFileOrString);
 *       else
 *         d_script=pcr_createScriptFromString(scriptFileOrString);
 *     }
 *     ~Script() {
 *       close();
 *     }
 *     public void close() {
 *       if (d_script!= IntPtr.Zero)
 *         pcr_destroyScript(d_script);
 *       d_script=IntPtr.Zero;
 *     }
 * 
 *     public String xmlReflection() {
 *       String s=pcr_ScriptXMLReflection(d_script);
 *       if (pcr_ScriptError(d_script)!=0)
 *         throw new Exception(d_script);
 *       close();
 *       return s;
 *     }
 *   }
 * 
 * public class ScriptFile : Script {
 *   public ScriptFile(String fileName) :
 *     base(fileName,true)
 *   {
 *   }
 *  }
 * public class ScriptString : Script {
 *   public ScriptString(String script String) :
 *     base(scriptString,false)
 *   {
 *   }
 * }
 */

 }
