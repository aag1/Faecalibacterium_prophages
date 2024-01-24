use strict;
use warnings;
use Bio::SeqIO;
use Bio::Annotation::Collection;



my $gbff_file = $ARGV[0];

my $seqio = Bio::SeqIO -> new(-file => $gbff_file, -format => "genbank" );



print "genome_id\tgenome_desc\tgenome_length\ttaxonomy\n";

while (my $entry = $seqio -> next_seq) {

    my $genome_id = $entry -> accession_number . "." . $entry -> version;

    my $genome_desc = $entry -> desc;

    my $genome_length = $entry -> length;

    my $taxonomy = join ";", reverse($entry -> species -> classification);

    print $genome_id . "\t" . $genome_desc . "\t" . $genome_length . "\t" . $taxonomy . "\n";

}
