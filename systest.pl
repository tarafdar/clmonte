#!/usr/bin/perl
system("./CLMonte");
system("./CLMonte_newSpin");
#test 123
my $file = 'outp.txt';
open my $info, $file or die "error no $file";

while (my $orig_line = <$info>) {
    my @orig_spin_array = split(' ', $orig_line);
    push(@orig_spin, $orig_spin_array[1]);
}



my $file = 'outp_newSpin.txt';
open my $info, $file or die "error no $file";

while (my $new_line = <$info>) {
    my @new_spin_array = split(' ', $new_line);
    push(@new_spin, $new_spin_array[1]);
}

close ($info);
$i=0;
foreach(@orig_spin){
    push(@sub_values, ($_ - $new_spin[$i]));     ;
    $i=$i + 1;
}

$i=0;
foreach(@sub_values){
    if($orig_spin[$i] != 0) 
    {
        push(@error_percentages, ($_/$orig_spin[$i]) * 100); 
    }
    else
    {
        push(@error_percentages, ($_/0.00001) * 100); 
    }
    $i=$i + 1;
}

open (OUTFILE, '>errorpcts.out');
$i=0;
foreach(@error_percentages){
    print OUTFILE "$i $_\n";
    $i=$i + 1;
}

close(OUTFILE);

open (OUTFILE, '>errorabs.out');
$i=0;
foreach(@sub_values){
    print OUTFILE "$i $_\n";
    $i=$i + 1;
}

close(OUTFILE);
