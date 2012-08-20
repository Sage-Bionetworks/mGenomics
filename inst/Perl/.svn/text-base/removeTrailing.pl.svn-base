my $f = shift;
my $newf = shift;

open F, $f;
$c = 0;
open O, ">".$newf;
while(<F>){
  $c++;
  chomp;
  my @a = split /\t/;
  my $i = scalar @a;
  my $pos = scalar @a - 1;
  while($i > 0){ 
    if($a[$i] !~ /[\w\d]/){ 
    #$a[$i] == "" or $a[$i] eq "" or $a[$i] != 0){ 
      $pos--;
    }else{
    #print $a[$i],"a$pos\t", join("\t", @a[0 .. $pos]), "\n";
#      print O $c, "\t", $pos, "\t", $a[10], "\t", $a[$pos], "\n";
      $pos+=1;
      my $string = join("\t", @a[0..$pos]). "\n";
      $string =~ s/["'#]//g;
      $string =~ s/\t\t/\tNA\t/g;
      print O $string;
      last;
    }
    $i--;
  }
}

close O;
