use strict;
use warnings;
use Bio::SeqIO;
use Bio::Annotation::Collection;
use Data::Dumper;



my $gpff_file = $ARGV[0];

my $seqio = Bio::SeqIO -> new(-file => $gpff_file, -format => "genbank" );



print "genome_id\tprotein_id\n";

while (my $entry = $seqio -> next_seq) {

    my $protein_id = $entry -> accession_number . "." . $entry -> version;


    my $genome_id = "";

    my @dblinks = $entry -> annotation -> get_Annotations("dblink");

    foreach my $x (@dblinks) {

        if ($x->{"database"} eq "REFSEQ") { $genome_id = $x->{"primary_id"} . "." . $x->{"version"} }

    }


    print $genome_id . "\t" . $protein_id . "\n";

}
