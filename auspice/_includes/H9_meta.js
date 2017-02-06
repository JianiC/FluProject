var	vaccineChoice = {};
vaccineChoice['A/Chicken/HongKong/G9/97'] = "1997-07-01";

var newChoice = {};
newChoice[''] = "";

var vaccineStrains = Object.keys(vaccineChoice);
var newStrains = Object.keys(newChoice);
var branch_labels= false;
var restrictTo = {"region": "all"};
//var restrictToh = {"host": "all"};

var genome_annotation = {'SP': [[1.2,1.2,1.2], [1,20,49]],
                         'HA1': [[1,1,1], [49,460,49+987]]};
var default_gene = 'HA1';

var structure = "5HMG.pdb"

var reference_viruses= {};
reference_viruses["A/Chicken/Beijing/1/94"] = '1994-07-01';
reference_viruses["A/chicken/Beijing/1/1994"] = '1994-07-01';
reference_viruses["A/Chicken/Korea/38349-p96323/96"] = '1996-07-01';
reference_viruses["A/Chicken/HongKong/G9/97"] = '1997-07-01';
reference_viruses["A/Duck/HongKong/Y280/97"] = '1997-07-01';
reference_viruses["A/Duck/HongKong/Y439/97"] = '1997-07-01';
reference_viruses["A/Quail/Hong_Kong/G1/97"] = '1997-07-01';

