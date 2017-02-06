var vaccineChoice = {};
vaccineChoice ['A/ruddyshelduck/Mongolia/1626/2010'] = "2010-09-11";

var newChoice = {};
newChoice[''] = "";

var vaccineStrains = Object.keys(vaccineChoice);
var newStrans = Object.keys(newChoice);
var branch_labels= false;
var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [1,20,52]],
                         'HA1':[[1,1,1], [52,460,52+981]],
                         'HA2':[[1.2,1.2,1.2], [52+981,1200,1701]]};
var default_gene = 'HA1';

var structure = "4LXV.pdb"

var reference_viruses = {};
