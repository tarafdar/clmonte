JedecChain;
	FileRevision(JESD32A);
	DefaultMfr(6E);

	P ActionCode(Ign)
		Device PartName(5SGXEA7H2) MfrSpec(OpMask(0));
	P ActionCode(Ign)
		Device PartName(EPM2210) MfrSpec(OpMask(0) SEC_Device(CFI_1GB) Child_OpMask(4 1 1 1 1) PFLPath("flash.pof"));

ChainEnd;

AlteraBegin;
	ChainType(JTAG);
AlteraEnd;
