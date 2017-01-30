
# TODO in bugzilla csharp topic #111

# just in case
CATEGORIES=notDefault
DLL_TARGET=PCRasterDotNetApi

# from devbin!
CSC=csc
REGASM="c:/WINDOWS/Microsoft.NET/Framework/v1.1.4322/regasm.exe"

first: rts


NUNITPATH=C:/projects/win32/NUnit-2.4.6-net-2.0/bin
NUNITlib=-lib:${NUNITPATH}/

NUNITrefs= -r:${NUNITPATH}/NUnit.Framework.dll
NUNIT=$(NUNITlib) $(NUNITrefs)


BUILD_DLL= -target:library -out:$(DLL_TARGET).dll
CSC_FLAGS= -unsafe $(BUILD_DLL) -define:DEBUG $(NUNIT)

PCRaster.cs: ../PCRaster.xsd $(XSD)/sources/xmlschema/commonTypes.xsd
	cp ../PCRaster.xsd $(XSD)/sources/xmlschema/commonTypes.xsd .
	xsd $< -c -n:pcrxml


CS_SRC= Script.cs PCRaster.cs XsdValidation.cs XsdValidationTest.cs
CS_SRC+= ScriptTest.cs

clean::
	rm -rf bin obj $(DLL_TARGET).exe $(DLL_TARGET).pdb
	rm -f *.xml *.xsd

all: $(CS_SRC)
	@$(CSC) $(CSC_FLAGS) $(CS_SRC)
	@cp PCRasterDotNetApi.dll $(QTDIR)/bin

# WILL break OpenMI test suite!
rts: all
	@cp sources/pcraster_model_engine/PCRaster.xsd .
	@cp ../testdata/apiExamples/*.xml .
	@${NUNITPATH}/nunit-console $(DLL_TARGET).dll
