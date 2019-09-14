using System;
using System.IO;
using System.Xml;
using System.Xml.Schema;

namespace pcr
{

  /// <summary>
  ///  support class to validate XML used in PCRaster API
  /// </summary>
  public class XsdValidation
  {
    bool          d_throwErrors;
    // int           d_nrMessages=0;
    String        d_filename;
    XmlTextReader d_xtr;


    public XsdValidation() {
      d_throwErrors=true;
    }
    public XsdValidation(bool throwErrors)
    {
      d_throwErrors=throwErrors;
    }

    private bool validateCommon() {
      bool success=true;
      //d_nrMessages=0;
      try {
      XmlValidatingReader vr = new XmlValidatingReader(d_xtr);

      vr.ValidationType = ValidationType.Schema;
      vr.ValidationEventHandler += new ValidationEventHandler (validationHandler);

      // consume
      while(vr.Read())
        ;
      } catch(Exception e) {
        success=false;
        if (d_throwErrors)
          throw e;
      } finally {
        d_xtr.Close();
      }
      return success;
    }

    public bool validateFile(String filename)
    {
      d_filename=filename;
      d_xtr = new XmlTextReader(d_filename);
      return validateCommon();
    }

    public bool validateString(String xmlString)
    {
      d_filename="?";

      StringReader  sr = new StringReader(xmlString);
      d_xtr = new XmlTextReader(sr);
      return validateCommon();
    }


    public void validationHandler(object sender, ValidationEventArgs args)
      {
        XmlSchemaException xse=args.Exception;
        // TODO dunno how xse encode the filename
        throw new Exception(String.Format(
         "{0}:{1}:{2}:{3}", d_filename,
           xse.LineNumber,xse.LinePosition, args.Message));
        // d_nrMessages +=1;
        // if (d_nrMessages == 3)
        // throw new Exception("Stopped at 3");
      }

   }
}
