// 2 color	["#5097BA", "#DF4327"]

var colors = [
	[],
	["#8EBC66"],
	["#4D92BF", "#E4662E"],
	["#4B8FC1", "#AABD52", "#E3612D"],
	["#4A8BC3", "#82BA72", "#CFB541", "#E25B2C"],
	["#4988C5", "#6EB389", "#AABD52", "#DEA73C", "#E2562B"],
	["#4785C7", "#64AD99", "#90BC65", "#C3BA46", "#E39A39", "#E1512A"],
	["#4682C9", "#5CA7A4", "#7FB975", "#AABD52", "#D2B340", "#E68F36", "#E04C29"],
	["#457FCB", "#57A1AD", "#73B584", "#96BD5F", "#BDBB49", "#DBAC3D", "#E68334", "#DF4628"],
	["#447BCD", "#539CB4", "#6AB090", "#88BB6C", "#AABD52", "#CBB842", "#E0A23A", "#E67A32", "#DF4127"],
	["#4377CD", "#5097BA", "#63AC9A", "#7CB879", "#9ABE5C", "#B9BC4A", "#D4B13F", "#E49938", "#E67030", "#DE3C26"],
	["#4273CE", "#4D93BE", "#5DA8A3", "#73B584", "#8DBC68", "#AABD52", "#C6B945", "#DBAC3D", "#E69036", "#E4672E", "#DD3725"],
	["#426FCE", "#4B8DC2", "#59A3AA", "#6BB18D", "#82BA71", "#9CBE5B", "#B7BD4B", "#CFB541", "#DFA43B", "#E68735", "#E35E2D", "#DD3124"]
];
var genotypeColors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"];

var epitopeColorScale = d3.scale.linear().clamp([true])
	.domain(epiColorDomain)
	.range(colors[10]);

var nonepitopeColorScale = d3.scale.linear().clamp([true])
	.domain(nonEpiColorDomain)
	.range(colors[10]);

var receptorBindingColorScale = d3.scale.linear().clamp([true])
	.domain(rbsColorDomain)
	.range(colors[4]);

var lbiColorScale = d3.scale.linear()
	.domain([0.0, 0.02, 0.04, 0.07, 0.1, 0.2, 0.4, 0.7, 0.9, 1.0])
	.range(colors[10]);

var dfreqColorScale = d3.scale.linear()
	.domain(dfreqColorDomain)
	.range(colors[10]);

var HIColorScale = d3.scale.linear()
	.domain(HIColorDomain)
	.range(colors[10]);

var cHIColorScale = d3.scale.linear()
	.domain(HIColorDomain)
	.range(colors[10]);

var dHIColorScale = d3.scale.linear().clamp([true])
	.domain(genericDomain.map(function (d){return 1.5*d;}))
	.range(colors[10]);

var regionColorScale = d3.scale.ordinal()
	.domain(regions.map(function(d){return d[0];}))
	.range(regions.map(function(d){return d[1];}));

var dateColorScale = d3.scale.linear().clamp([true])
	.domain(dateColorDomain)
	.range(colors[10]);
	
/*var hostColorScale = d3.scale.ordinal()
	.domain(hosts.map(function(d){return d[0];}))
	.range(hosts.map(function(d){return d[1];}));*/
	
var regionColorScale = d3.scale.ordinal()
	.domain(regions.map(function(d){return d[0];}))
	.range(regions.map(function(d){return d[1];}));

/*var hostColorScale = d3.scale.ordinal()
	.domain(hosts.map(function(d){return d[0];}))
	.range(hosts.map(function(d){return d[1];}));*/

var fitnessColorScale = d3.scale.linear().clamp([true])
	.domain(fitnessColorDomain)
	.range(colors[10]);

// "ep", "ne" and "rb" need no adjustments
function adjust_coloring_by_date() {
	if (colorBy == "lbi") {
		calcLBI(rootNode, nodes, false);
		nodes.forEach(function (d) {
			d.coloring = d.LBI;
		});
	}
	else if (colorBy == "date") {
		nodes.forEach(function (d) {
			d.coloring = d.num_date;
		});
	}
}

function stateAtPosition(clade, gene, pos){
	if (typeof cladeToSeq[clade][gene][pos] == "undefined"){
		return cladeToSeq["root"][gene][pos];
	}else{
		return cladeToSeq[clade][gene][pos];
	}
}

function colorByTrait() {

	colorBy = document.getElementById("coloring").value;
	if (colorBy=="--"){
		document.getElementById("coloring").value = "ep";
		colorBy = document.getElementById("coloring").value;
	}
	console.log(colorBy);
	d3.selectAll('.serum')
		.style("visibility", serumVisibility);
	var vis = (colorBy=='HI_dist')?'block':'none';
	if (document.getElementById("HIcontrols") !== null) {
			document.getElementById("HIcontrols").style.display = vis;
	}

	if (colorBy == "ep") {
		colorScale = epitopeColorScale;
		nodes.map(function(d) { d.coloring = d.ep; });
	}
	else if (colorBy == "ne") {
		colorScale = nonepitopeColorScale;
		nodes.map(function(d) { d.coloring = d.ne; });
	}
	else if (colorBy == "rb") {
		colorScale = receptorBindingColorScale;
		nodes.map(function(d) { d.coloring = d.rb; });
	}
	else if (colorBy == "lbi") {
		colorScale = lbiColorScale;
		adjust_coloring_by_date();
	}
	else if (colorBy == "dfreq") {
		colorScale = dfreqColorScale;
		nodes.map(function(d) { d.coloring = d.dfreq;});
	}
	else if (colorBy == "region") {
		colorScale = regionColorScale;
		nodes.map(function(d) { d.coloring = d.region; });
	}
	else if (colorBy == "cHI") {
		colorScale = cHIColorScale;
		nodes.map(function(d) { d.coloring = d.cHI; });
	}
	else if (colorBy == "date") {
		colorScale = dateColorScale;
		nodes.map(function(d) { d.coloring = d.num_date; });
	}
	//else if (colorBy == "host") {
		//colorScale = hostColorScale;
		//nodes.map(function(d) { d.coloring = d.host; });
	//}
	else if (colorBy == "fitness") {
		colorScale = fitnessColorScale;
		nodes.map(function(d) { d.coloring = d.fitness; });
	}

	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);

	d3.selectAll(".tip")
		.attr("r", tipRadius)
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);

	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();
}

function tipStrokeColor(d) {
	var col = colorScale(d.coloring);
	return d3.rgb(col).toString();
}

function tipFillColor(d) {
	var col = colorScale(d.coloring);	;
	return d3.rgb(col).brighter([0.65]).toString();
}

function branchStrokeColor(d) {
	var col;
	if (colorBy == "region" || colorBy == "date") {
		col = "#AAA";
	}
	else {
		if (typeof d.target.coloring != "undefined"){
			col = colorScale(d.target.coloring);
		}else{
			col="#AAA";
		}
	}
	var modCol = d3.interpolateRgb(col, "#BBB")(0.6);
	return d3.rgb(modCol).toString();
}

//function branchStrokeColor(d) {
//	var col;
//	if (colorBy == "host" || colorBy == "date") {
//		col = "#AAA";
//	}
//	else {
//		if (typeof d.target.coloring != "undefined"){
//			col = colorScale(d.target.coloring);
//		}else{
//			col="#AAA";
//		}
//	}
//	var modCol = d3.interpolateRgb(col, "#BBB")(0.6);
//	return d3.rgb(modCol).toString();
//}

function contains(arr, obj) {
    for(var i=0; i<arr.length; i++) {
        if (arr[i] == obj) return true;
    }
}

function parse_gt_string(gt){
	mutations = [];
	gt.split(',').map( function (d) {
		var tmp = d.split(/[\s//]/); //FIXME: make more inclusive
		var region;
		var positions = [];
		for (var i=0; i<tmp.length; i++){
			if (contains(["EU","NA","AS","OC"], tmp[i])){
				region = tmp[i];
			}else{
				if (tmp[i].length>0) positions.push(tmp[i]);
			}
		}
		if (typeof region == "undefined") region="global";
		// sort if this is a multi mutation genotype
		if (positions.length>1){
			positions.sort(function (a,b){
				return parseInt(a.substring(0,a.length-1)) - parseInt(b.substring(0,b.length-1));
			});
		}
		mutations.push([region, positions.join('/')]);
	});
	return mutations;
};

function colorByGenotype() {
	var positions_string = document.getElementById("gt-color").value.split(',');
	var positions_list = []
	positions_string.map(function(d) {
		var pos_fields = d.split(':');
		var val, gene;
		if (pos_fields.length==1){
			val = parseInt(pos_fields[0])-1;
			gene=default_gene;
		}else if (pos_fields.length==2){
			val = parseInt(pos_fields[1])-1;
			gene=pos_fields[0].replace(' ','');
		}else{
			val = parseInt('NaN');
		}
		console.log('attempt genotype coloring: '+ [gene, val]);
		if ((!isNaN(val))&&(typeof cladeToSeq["root"][gene]!="undefined")) {
			if (val < cladeToSeq["root"][gene].length) {
				positions_list.push([gene, val]);
			}
		}
	});
	console.log(positions_list);
	if (positions_list.length > 0) {
		colorBy = "genotype";
		colorByGenotypePosition(positions_list);
	}
	else {
		d3.select("#coloring").each(colorByTrait);
		gt = parse_gt_string(freqdefault);
		make_gt_chart(gt);
		document.getElementById("gtspec").value = freqdefault;
	}
}


function colorByGenotypePosition (positions) {
	var gts = nodes.map(function (d) {
		var tmp = [];
		for (var i=0; i<positions.length; i++){
			tmp[tmp.length] = positions[i][0]+':'+(positions[i][1]+1)+stateAtPosition(d.clade, positions[i][0], positions[i][1]);
		}
		d.coloring = tmp.join('/');
		return d.coloring;});
	var unique_gts = d3.set(gts).values();
	var gt_counts = {};
	for (var i=0; i<unique_gts.length; i++){gt_counts[unique_gts[i]]=0;}
	gts.forEach(function (d) {gt_counts[d]+=1;});
	var filtered_gts = unique_gts.filter(function (d) {return gt_counts[d]>=10;});
	filtered_gts.sort(function (a,b){
		var res;
		if (gt_counts[a]>gt_counts[b]){ res=-1;}
		else if (gt_counts[a]<gt_counts[b]){ res=1;}
		else {res=0;}
		return res;});
	console.log("genotypes passed filtering:"+filtered_gts);
	colorScale = d3.scale.ordinal()
		.domain(filtered_gts)
		.range(genotypeColors);
	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);
	treeplot.selectAll(".tip")
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);
	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();

	if ((positions.length==1)&&(filtered_gts.length>1)){
		var tmp_gts=[];
		for (var ii=0; ii<filtered_gts.length; ii+=1){
			tmp_gts.push(["global", filtered_gts[ii]])
		}
		make_gt_chart(tmp_gts);
		document.getElementById("gtspec").value = tmp_gts.map( function (d) {return d[1];}).join(', ');
	}
}


function newFocus(){
	if (typeof(focusNode)=="undefined"){
		var ntiters = 0, ntmp;
		focusNode=sera[0];
		for (var i=0; i<sera.length; i++){
			ntmp = Object.keys(sera[i].mean_HI_titers).length;
			if (ntmp>ntiters){
				ntiters = ntmp;
				focusNode = sera[i];
			}
		}
	}
	// add checkboxes to include/exclude sera
	var seraDiv = document.getElementById("sera");
	var htmlStr = "";
	activeSera = {};
	console.log(focusNode);
	for (var serum in focusNode.potency_mut){
		var serumID = serum.split("/").join("");
		htmlStr+='<input type="checkbox" id="' + serumID + '" name="' + serum + '" checked="checked"> ' + serum +"<br>";
		activeSera[serum]=true;
	}
	seraDiv.innerHTML = htmlStr;
	console.log(seraDiv);
	for (var serum in focusNode.potency_mut){
		var serumID = serum.split("/").join("");
		d3.select("#"+serumID)
			.on("change", function(elem){
					for (var tmpserum in focusNode.potency_mut){
						var tmpserumID = tmpserum.split("/").join("");
						activeSera[tmpserum]=document.getElementById(tmpserumID).checked;
					}
					colorByHIDistance()});
		}

	colorByHIDistance();
}

function colorByHIDistance(){
	correctVirus = document.getElementById("virus").checked;
	correctPotency = document.getElementById("serum").checked;
	var HIchoices = document.getElementsByName("HImodel");
	console.log(activeSera);
	for(var i = 0; i < HIchoices.length; i++){
	    if(HIchoices[i].checked){
	        HImodel = HIchoices[i].value;
	    }
	}
	colorBy = 'HI_dist'

	treeplot.selectAll(".serum")
	.style("fill", function (d){if (d==focusNode) {return '#FF3300';} else {return '#555555';}})
		.style("font-size", function (d) {if (d==focusNode) {return "30px";} else {return "12px";}})
		.text(function (d) {if (d==focusNode) {return '\uf05b';} else {return serumSymbol;}});

	console.log("Using HI model: "+HImodel);
	console.log("Color by HI Distance from "+focusNode.strain);
	console.log("correcting for virus effect: "+correctVirus);
	console.log("correction for serum effect: "+correctPotency);

	calcHImeasured(focusNode, rootNode);
	calcHImutations(focusNode, rootNode);
	calcHItree(focusNode, rootNode);

	colorScale = HIColorScale;
	if (HImodel=='mutation'){
		nodes.map(function(d) { d.coloring = d.HI_dist_mut;});
	}else if (HImodel=='tree'){
		nodes.map(function(d) { d.coloring = d.HI_dist_tree;});
	}else{
		nodes.map(function(d) { d.coloring = d.HI_dist_meas;});
	}

	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);

	treeplot.selectAll(".tip")
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);

	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();
}


d3.select("#coloring")
	.style("cursor", "pointer")
	.on("change", colorByTrait);


var genotypeColoringEvent;
d3.select("#gt-color")
	.on("keyup", function(){
		if (typeof genotypeColoringEvent != "undefined"){clearTimeout(genotypeColoringEvent);}
		genotypeColoringEvent = setTimeout(colorByGenotype, 200);
	});



