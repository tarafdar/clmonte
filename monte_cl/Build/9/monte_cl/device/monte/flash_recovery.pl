#This is a recovery script file for the Nalla_pcie BSPs

my $aocxfile = $ARGV[0];


if ( $aocxfile =~ m/\.aocx$/i ) {
  -f $aocxfile or print "Error: can't find $aocxfile file needed for flashing." and exit 1;
} else {
  print "Error: Currently only .aocx files are accepted.\n";
  print "       Please use the .aocx from the original compiled project or the base.aocx included with the target BSP\n";
  exit 1;
}

my $binfile = "fpga_temp.bin";
my $soffile = "fpga_temp.sof";
my $pci_id_file = "pci_id_temp.txt";

system ("aocl binedit $aocxfile get .acl.fpga.bin $binfile");
$? == 0  or die "Error: aocl binedit failed getting fpga_temp.bin";

system ("aocl binedit $binfile get .acl.sof $soffile");
$? == 0  or die "Error: aocl binedit failed getting fpga_temp.sof";

$pci_id = `aocl binedit $binfile print .acl.pcie.dev_id`;

if ( $pci_id =~ /A700/  || $pci_id =~ /D500/ ) {
    print "PCIe385n AOCX file\n";

my $cof = <<END;
<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<cof>
	<eprom_name>CFI_1GB</eprom_name>
	<output_filename>flash.pof</output_filename>
	<n_pages>2</n_pages>
	<width>1</width>
	<mode>12</mode>
	<sof_data>
		<start_address>02140000</start_address>
		<user_name>Page_2</user_name>
		<page_flags>4</page_flags>
		<bit0>
			<sof_filename>$soffile</sof_filename>
		</bit0>
	</sof_data>
	<sof_data>
		<start_address>00040000</start_address>
		<user_name>Page_1</user_name>
		<page_flags>2</page_flags>
		<bit0>
			<sof_filename>$soffile</sof_filename>
		</bit0>
	</sof_data>
	<version>5</version>
	<create_cvp_file>0</create_cvp_file>
	<options>
		<map_file>1</map_file>
		<option_start_address>0</option_start_address>
		<dynamic_compression>0</dynamic_compression>
	</options>
</cof>
END

my $cdf = <<CDFEND;
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
CDFEND

open COFFILE, ">flash.cof";
print COFFILE $cof;
close COFFILE;

open CDFFILE, ">flash.cdf";
print CDFFILE $cdf;
close CDFFILE;

system ("quartus_cpf --convert flash.cof");
$? == 0  or die "Error: quartus_cpf failed";

print "The following will attempt to setup the Nallatech USB-Blaster-II Breakout & Debug board if connected\n";
print "Note if this is not present please ignore any error messages regarding setting this cable up\n";
print "Setting FPGA Configurator Jtag clock rate to: \n";
system ('jtagconfig','--setparam',"Nallatech PCI3-395",'JtagClock','16M');
system ('jtagconfig','--getparam',"Nallatech PCI3-395",'JtagClock');
print "Finished attempting to setup Nallatech USB-Blaster-II Breakout & Debug board\n";

system ("quartus_pgm -c 1 flash.cdf");
$? == 0  or die "Error: quartus_cpf failed";

} elsif( $pci_id =~ /D800/ || $pci_id =~ /AB00/  ) {
    print "PCIe395 AOCX file\n";

my $cof = <<END;
<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<cof>
	<eprom_name>CFI_512Mb</eprom_name>
	<output_filename>flash.pof</output_filename>
	<n_pages>1</n_pages>
	<width>1</width>
	<mode>12</mode>
	<sof_data>
		<start_address>10000</start_address>
		<user_name>Page_1</user_name>
		<page_flags>2</page_flags>
		<bit0>
			<sof_filename>$soffile</sof_filename>
		</bit0>
	</sof_data>
	<version>5</version>
	<create_cvp_file>0</create_cvp_file>
	<options>
		<map_file>1</map_file>
		<option_start_address>0</option_start_address>
		<dynamic_compression>0</dynamic_compression>
	</options>
</cof>
END

my $cdf = <<CDFEND;
JedecChain;
	FileRevision(JESD32A);
	DefaultMfr(6E);

	P ActionCode(Ign)
		Device PartName(EPM2210) MfrSpec(OpMask(0) SEC_Device(QSPI_512MB) Child_OpMask(3 1 1 1) PFLPath("flash.pof"));
	P ActionCode(Ign)
		Device PartName(5SGSMD8N2) MfrSpec(OpMask(0));

ChainEnd;

AlteraBegin;
	ChainType(JTAG);
AlteraEnd;


CDFEND

open COFFILE, ">flash.cof";
print COFFILE $cof;
close COFFILE;

open CDFFILE, ">flash.cdf";
print CDFFILE $cdf;
close CDFFILE;

system ("quartus_cpf --convert flash.cof");
$? == 0  or die "Error: quartus_cpf failed";

system ("quartus_pgm -c 1 flash.cdf");
$? == 0  or die "Error: quartus_cpf failed";

} else {
    print "Unknown AOCX file\n"
    }

exit
