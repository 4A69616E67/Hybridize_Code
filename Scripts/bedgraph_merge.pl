###########bedgraph format as input file
###########merged the adjacent bins and calculate the average RPKM value

$in_file = $ARGV[0];
open(IN,$in_file) or die "$!\n";
my($chr,$start,$end,$score,$len);
$line=<IN>;
chomp $line;
@strs = split(/\s+/,$line);
$chr = $strs[0];
$start = $strs[1];
$end = $strs[2];
$len = $end-$start;
$score = $strs[3]/10*$len;
while($line=<IN>){
	chomp $line;
	@strs = split(/\s+/,$line);
	###判断是毗邻
	if($chr eq $strs[0] && $end == $strs[1]){
		$end =$strs[2];
		$len = $end - $start;
		$score += $strs[3]/10*($strs[2]-$strs[1]);
	}else{
		print $chr."\t".$start."\t".$end."\t".($score/$len*10)."\n";
		$chr = $strs[0];
		$start = $strs[1];
		$end = $strs[2];
		$len = $end-$start;
		$score = $strs[3]/10*$len;
	}
}
print $chr."\t".$start."\t".$end."\t".($score/$len*10)."\n";