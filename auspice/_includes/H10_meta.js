var	vaccineChoice = {};
vaccineChoice['A/mallard/Maryland/13OS3318/2014'] = "2014-06-24";

var	newChoice = {};
newChoice['A/Unknown/Unknown/Batch2-1_002_01102017_4_H10N6'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_007_01102017_4_H10N4'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_008_01102017_4_H10N7'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_010_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_013_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_021_01102017_4_H10N4'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_027_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_028_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_030_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_031_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_033_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_034_01102017_4_H10N7'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_036_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_037_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_039_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_042_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_043_01102017_4_H10N4'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_044_01102017_4_H10N5'] = "2017-01-10";
newChoice['A/Unknown/Unknown/Batch2-1_048_01102017_4_H10N5'] = "2017-01-10";

var vaccineStrains = Object.keys(vaccineChoice);
var newStrains = Object.keys(newChoice);
var branch_labels= false;

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [12,20,57]],
                         'HA1':[[1,1,1], [57,460,57+1038]],
                         'HA2':[[1.2,1.2,1.2], [57+1038,1200,1769]]};
var default_gene = 'HA1';

var structure = "4FQM.pdb";

var reference_viruses = {};