
sub commit_utils($)
{
	my $build_dir = $_[0];
	system("cp -r $ENV{'MONTE_ROOT'}/utils $build_dir/.");
	return 1;
}

sub commit_test($)
{
	my $build_dir = $_[0];
	system("cp -r $ENV{'MONTE_ROOT'}/test $build_dir/.");
	return 1;
}

sub commit_makefile($)
{
	my $build_dir = $_[0];
	system("cp $ENV{'MONTE_ROOT'}/Makefile $build_dir/.");
	return 1;
}

sub commit_bin($)
{
	my $build_dir = $_[0];
	system("cp -r $ENV{'MONTE_ROOT'}/bin $build_dir/.");
	return 1;
}

sub commit_source($)
{
	my $build_dir = $_[0];
	system("cp -r $ENV{'MONTE_ROOT'}/monte_cl $build_dir/.");
	return 1;
}

sub make_build_dir($)
{
	my $build_num = $_[0];
	my $build_dir = "$ENV{'MONTE_ROOT'}/Build/$build_num";
	mkdir $build_dir;
	return $build_dir;
}

sub update_build_number($)
{
	my $build_num = $_[0];
	my $fh;
	open($fh, "+>", "$ENV{'MONTE_ROOT'}/meta/meta.txt");
	print $fh "build_num=$build_num";
	close $fh;
}

sub get_next_build_number()
{
	my $build_num = -1;

	my $fh;
	open($fh, "+<", "$ENV{'MONTE_ROOT'}/meta/meta.txt");
	my @lines = <$fh>;
	foreach my $line (@lines) {
		my @parts = split(/=/, $line, 2);
		my $key = $parts[0];
		my $value = $parts[1];

		print "Key: $key\n";
		print "Value: $value\n";
		if($key eq "build_num") {
			$build_num = $value;
		}
	}

	close $fh;
	return $build_num + 1;
}

if(!defined($ENV{'MONTE_ROOT'})) {
	print "MONTE_ROOT not set!\n";
	return 0;
}

my $build_num = get_next_build_number();
$build_num = $build_num;
my $build_dir = make_build_dir($build_num);
commit_source($build_dir);
commit_bin($build_dir);
commit_test($build_dir);
commit_utils($build_dir);
commit_makefile($build_dir);
update_build_number($build_num);

print "Commit $build_dir Complete!\n";

1;
