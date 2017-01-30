using System;
using NUnit.Framework;


namespace pcr.Test {


/// <summary>
///  Test class to demonstrates the C-API
///  OO-functionality like ScriptFile or ScriptString are not used
///  so the examples are easy to port to other languages.
/// </summary>
[TestFixture]
public class XsdValidationTest {

  static String document(String docElName,String docElAttrs,String innerXml) {
   return String.Format(
    @"<{0}
      xmlns='http://www.pcraster.nl/pcrxml'
      xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'
      xsi:schemaLocation='http://www.pcraster.nl/pcrxml PCRaster.xsd'
      {1}
      >
      {2}
      </{3}>",docElName,docElAttrs,innerXml,docElName);
  }
  static String definitionDocument(String attrs,String innerXml) {
    return document("definition",attrs,innerXml);
  }
  // expect no errors, throw if we have errors
  static XsdValidation validValidator = new XsdValidation(true);
  static XsdValidation badValidator   = new XsdValidation(false);

  [Test] public void TestDefinition()
  {
   Assert.IsTrue(validValidator.validateString(definitionDocument("name='a'","")));

   String xml = definitionDocument("name='a'","<field/>");
   Assert.IsTrue(validValidator.validateString(xml));
  }
  [Test] public void TestNotValid()
  {
   Assert.IsFalse(badValidator.validateFile("notValid.xml"));
  }
  [Test] public void TestApiExample1()
  {
   Assert.IsTrue(validValidator.validateFile("lookup.xml"));
  }
  [Test] public void TestXMLBindig()
  {
    pcrxml.Script s= pcr.API.xmlFileToXMLScript("lookup.xml");
    Assert.IsNull(s.definition[0].Item);
    Assert.AreEqual("inLineTable",s.definition[2].name);

    // optional element
    Assert.IsNotNull(s.definition[0].scriptInput);
    Assert.IsTrue(s.definition[0].scriptInput != null);
    Assert.IsNull(s.definition[0].scriptOutput);
    Assert.IsTrue(s.definition[0].scriptOutput == null);

    // choice element
    Assert.IsNotNull(s.definition[2].Item);
    Assert.IsTrue(False); // bugzilla #111
    // Assert.IsTrue(s.definition[2].Item is pcrxml.Relation);
    // Assert.IsFalse(s.definition[2].Item is pcrxml.Field);
  }
}

} // namespace pcr.XsdValidationTest
