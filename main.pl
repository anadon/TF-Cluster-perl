use LWP::Simple;
use File::Copy;
use Cwd;
# use Config::Simple;
# use File::Basename;

my $username = shift;
# my $directory = shift;
my $basename = shift;

my $uploadFilesPath = shift;


my $cutoff = shift;
my $theta1 = shift;
my $theta2 = shift;
my $theta3 = shift;
my $corr_name = shift;
my $tf_list = shift;
my $expressionData = shift;
my $decom = shift;
my $SSGA_p = shift;
my $SSGA_t = shift;
my $MSGA_p = shift;
my $MSGA_q = shift;
my $MSGA_t = shift;


my $start = time;


if ($corr_name == 3 ) {
	# Weighted Rank Correlation
	`perl /var/www/html/cluster/perl/SCCM_v1.1/scripts/SCCM_pipe.pl ./SCCM_pipe.cfg`;
}
else  {
	# run the c++ program
	my $filename = 'report.txt';
	if ($corr_name == 1) {
		$corr_name = "spearman";
	} elsif ($corr_name == 2) {
		$corr_name = "pearson";
	} elsif ($corr_name == 4) {
		$corr_name = "kendall";
	} else {
		$corr_name = "maximum-information-coefficient";
	}
	open(my $fh0, '>', $filename) or die "Could not open file '$filename' $!";
	my $tf_listWholePath = $uploadFilesPath . "/" . $tf_list;
	my $expressionDataWholePath = $uploadFilesPath . "/" . $expressionData;
	# print $fh0 `/var/www/html/cluster/cplusplus/tf-cluster/TF-cluster/tf-cluster -1 1.5 -2 1.2 -3 0.8 -c spearman -e ${tf_listWholePath} -t ${expressionDataWholePath} -k 100`;
	print $fh0 `/var/www/html/cluster/cplusplus/tf-cluster/tf-cluster -1 ${theta1} -2 ${theta2} -3 ${theta3} -c ${corr_name} -e ${expressionDataWholePath} -t ${tf_listWholePath} -k ${cutoff}`;

	close $fh0;
	# Separate the report.txt to two files. One is SCCM matrix. Another is cluster output.
	# These two files separate with '---###???'
	`csplit --digits=2 --quiet --prefix=outfile report.txt "/---###???/+1" "{*}"`;

	# get file names
	my $outfileSCCM = $tf_list . "_" . $expressionDataName . "_coexp.matrix.txt";


	# remove the separate characters of '---###???' at the last 5 lines of the first file
	qx(tac outfile00 | sed '1,4 d' | tac > $outfileSCCM);

	if ($decom eq "SSGA") {
		`rm outfile01`;
		qx(perl /var/www/html/cluster/perl/SCCM_v1.1/scripts/SSGA.pl -m $outfileSCCM -p $SSGA_p -t $SSGA_t >>standout);
	} elsif ($decom eq "MSGA") {
		`rm outfile01`;
		qx(perl /var/www/html/cluster/perl/SCCM_v1.1/scripts/MSGA.pl -m $outfileSCCM -p $MSGA_p -q $MSGA_q -t $MSGA_t>>standout);
	} else {
		my $outfileCluster = "Cluster_" . $tf_list . "_" . $expressionDataName . "_coexp.matrix.txt";
		`cp outfile01 ${outfileCluster}`;
	}
}

 my $duration = time - $start;

 my $run_time = sprintf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($duration / 3600), int(($duration % 3600) / 60), int($duration % 60));

my $curPath = cwd();

print $curPath;
@files = <*>;
foreach $file (@files) {
  if (index($file, "matrix.txt") != -1) {
	move("./$file", "../running_result/$file");
  }
}


chdir('../running_result');
my $zip_output = 'allOutResult.zip';



system('zip -r '.$zip_output.' ./');
#Send email to the user, attach the result

my $attachmentpath = "User/$username/$basename/running_result/allOutResult.zip";

my $url = "http://sys.bio.mtu.edu/cluster/sendemail.php?usr=$username&file=$attachmentpath&run_time=$run_time";

open(my $fh, '>', 'report.txt');
print $fh "$run_time\n";
print $fh "done\n";
close $fh;
# print $url;
# my $url = "http://sys.bio.mtu.edu/cluster/sendemail.php?usr=$username&basename=$basename";
getprint($url);



