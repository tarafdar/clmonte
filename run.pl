$input_file = "input.txt";
$CLMonte_outFile = "output_newSpin.txt";
$plot_pic = "newPlot.png";

open my $info, $input_file or die "Could not open $file ";

#input parameters = (G, MUS_MAX, V, COS_CRIT, N);
@input_params = ("0.9" , "90.0" , "0.0214" , "0.6999", "1.4"); 
while(my $line = <$info>){
    my @values = split(' ', $line);
    if($values[0] eq 'G'){
        $input_params[0] = $values[1];
     }
    elsif($values[0] eq 'MUS_MAX'){
        $input_params[1] = $values[1];
     }
    elsif($values[0] eq 'V') {
        $input_params[2] = $values[1];
     }
    elsif($values[0] eq 'COS_CRIT'){
        $input_params[3] = $values[1];
     }
    elsif($values[0] eq 'N'){
        $input_params[4] = $values[1];
     }
}

print "input parameters\nG = $input_params[0] \nMUS_MAX = $input_params[1] \nV = $input_params[2] \nCOS_CRIT = $input_params[3] \nN = $input_params[4]\n";
system("./CLMonte_newSpin $input_params[0] $input_params[1] $input_params[2] $input_params[3] $input_params[4] $CLMonte_outFile");
system("gnuplot -e \"filename = \'$CLMonte_outFile\';plot_out=\'$plot_pic\' ;\" plot.pt"); 
