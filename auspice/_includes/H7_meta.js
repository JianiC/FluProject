var	vaccineChoice = {};
vaccineChoice['A/turkey/Italy/3889/99'] = "1999-07-01";
vaccineChoice['A/mallard/Netherlands/12/00'] = "2000-07-01";
vaccineChoice['A/human/New-York/107/2003'] = "2003-07-01";
vaccineChoice['A/human/Shanghai/2/2013'] = "2013-03-05";
vaccineChoice['A/human/Anhui/1/2013'] = "2013-03-20";

var newChoice = {};
newChoice[''] = "";

var vaccineStrains = Object.keys(vaccineChoice);
var newStrains = Object.keys(newChoice);
var branch_labels= false;
var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [12,20,57]],
                         'HA1':[[1,1,1], [57,460,57+1038]],
                         'HA2':[[1.2,1.2,1.2], [57+1038,1200,1769]]};
var default_gene = 'HA1';

var structure = "4M40.pdb";

var reference_viruses = {};

