
function getQueryVariable(variable) {
	var query = window.location.search.substring(1);
	var vars = query.split("&");
	for (var i=0;i<vars.length;i++) {
		var pair = vars[i].split("=");
		if (pair[0] == variable) {
			return pair[1];
		}
	}
}

var pi = Math.PI;

function overlaps(s1, e1, s2, e2) {
    if (s1 < s2 & e1 < e2 & e1 > s2) {
        return 1;
    } else if (s2 < s1 & e2 > s1 & e2 < e1) {
        return 1;
    } else if (s2 > s1 & s2 < e1 & e2 < e1) {
        // return 1;
    } else if (s2 < s1 & e2 > e1) {
        return 1;
    } else if (s1 < s2 & e1 > e2) {
        return 1;
    } else if (s1 == s2 | s1 == e2 | e1 == s1 | e1 == e2) {
        return 1;
    }
    return 0;
}


function adjustBump(annot, offset) {
    var recurse = 0;
    for (var i = 0; i < annot.length; i++) {
        for (var j = 0; j < annot.length; j++) {
            if (i != j) {
                var g1 = annot[i];
                var g2 = annot[j];
                if (g1.bumpLevel == g2.bumpLevel && overlaps(g1.start - offset, g1.end + offset, g2.start - offset, g2.end + offset)) {
                    annot[i].bumpLevel++;
                    recurse = 1;
                }
            }
        }
    }
    if (recurse) {
        adjustBump(annot, offset);
    }
}


function computeCartesian(r, coord, totalBP) {
    var arcAvail = 360 - 10;
    var ratio = coord / totalBP;
    var theta = (((coord / totalBP) * arcAvail) * (pi / 180)) + (5 * (pi / 180));
    if (theta <= pi / 2) {
        return ({
            x: r * Math.sin(theta),
            y: r * Math.cos(theta) * -1
        });
    } else if (theta > pi / 2 && theta <= pi) {
        return ({
            x: r * Math.sin(theta),
            y: r * Math.cos(theta) * -1
        });
    } else if (theta > pi && theta <= (3 * pi) / 2) {
        return ({
            x: r * Math.sin(theta),
            y: r * Math.cos(theta) * -1
        });
    } else if (theta > (3 * pi) / 2 && theta <= 2 * pi) {
        return ({
            x: r * Math.sin(theta),
            y: r * Math.cos(theta) * -1
        });
    } else {
        theta = (arcAvail * (pi / 180)) + (5 * (pi / 180))
        return ({
            x: r * Math.sin(theta),
            y: r * Math.cos(theta) * -1
        });
    }


}


function computePath(start, end, r, totalBP, diameter) {
    // creates some d magic to connect paths
    // <path class="SamplePath" d="M100,200 C100,100 250,100 250,200
    //                                 S400,300 400,200" />
    startcoords = computeCartesian(r, start, totalBP);
    endcoords = computeCartesian(r, end, totalBP);
    //harcoded !!!!!!!!
    startcontrol = computeCartesian(r - (diameter * 0.1), start, totalBP);
    endcontrol = computeCartesian(r - (diameter * 0.1), end, totalBP);
    return ("M" + startcoords.x + "," + startcoords.y +
        " C" + startcontrol.x + "," + startcontrol.y + "," + endcontrol.x + "," + endcontrol.y + " " + endcoords.x + "," + endcoords.y);
}



function computeStrandPath(start, end, r, totalBP, flag) {
    startcoords = computeCartesian(r, start, totalBP);
    endcoords = computeCartesian(r, end, totalBP);
    //var flag = "0,1";
    if ((end - start) /totalBP > 0.5){
    //flag = "1,1";
    //flag = "0,0";
    }
    return ("M" + startcoords.x + "," + startcoords.y +
        " A" + r + "," + r + " 0 " + flag + " " + endcoords.x + "," + endcoords.y);
}

function computeArcPath(start, end, r1, r2, totalBP) {
    startcoords1 = computeCartesian(r1, start, totalBP);
    endcoords1 = computeCartesian(r1, end, totalBP);
    startcoords2 = computeCartesian(r2, start, totalBP);
    endcoords2 = computeCartesian(r2, end, totalBP);
    var flag1 = "0,1";
    if ((end - start) /totalBP > 0.5){
        	flag1 = "1,1";
        }
        var flag2 = "0,0";
        if ((end - start) /totalBP > 0.5) {
        flag2 = "0,1";
    }
    return ("M" + startcoords1.x + "," + startcoords1.y +
        " A" + r1 + "," + r1 + " 0 " + flag1 + " " + endcoords1.x + "," + endcoords1.y +
        " L" + endcoords2.x + "," + endcoords2.y +
        " A" + r2 + "," + r2 + " 0 " + flag2 + " " + startcoords2.x + "," + startcoords2.y +
        " z");
}

function computePointPath(start, end, score, minscore, maxscore, r, totalBP, diameter) {
    //var maxscore = 20;
    var adjMaxscore = maxscore - minscore;
    var adjScore = score - minscore;
    var trackwidth = diameter * 0.05;
    var radius = r + ((parseFloat(adjScore) / adjMaxscore) * trackwidth);
    var startcoords = computeCartesian(radius, start, totalBP);
    return "translate(" + (startcoords.x + (diameter * 0.5)) + "," + (startcoords.y + (diameter * 0.5)) + ")";
}

function numberWithCommas(x) {
    return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}




function renderHic(gene, tissue, diameter, breadcrumb) {
    resetPage(gene, tissue, breadcrumb)
    var trans = "translate(" + diameter * 0.5 + "," + diameter * 0.5 + ")";
    var maxscore = 0;
    d3.json("/CHIC_DEMO2/cgi-bin/prototype.pl?gene=" + gene + '&tissue=' + tissue, function (error, json) {
        if (error) return console.warn(error);
        data = json;
        //console.log(data);
        var genes = data.genes;
        var snps = data.snps;
        var meta = data.meta;
        //console.log(data.region);
        // compute max score
        //var maxscore = 0;
        for (var i = 0; i < snps.length; i++) {
            if (snps[i].score > maxscore) {
                maxscore = parseFloat(snps[i].score);
            }
        }
        //console.log(snps);
        d3.select("#message").remove();
        var hics = data.hic;
        if (hics.length == 0) {
            div = d3.select("#svg-container")
                .append("div")
                .html("<h1>No interactions found</h1>")
                .attr("id", "message")
                .style("width", "100%")
                .style("text-align", "center")
                .style("padding-top", "200px");
            return;
        }
        for (var i = 0; i < hics.length; i++) {
            hics[i].id = i + 1;
        }
        // set this to make genes that are close but not overlapping bump
        var offset = 0;
        adjustBump(genes, offset);
        var bt = {};
        for (var g in genes) {
            bt[genes[g].color] = 1;

        }
        bt['hilight'] = 1;
        var bt_scale = d3.scale.category10().domain(Object.keys(bt));
        var totalBP = data.meta.rend - data.meta.rstart;
        var vis = d3.select("#svg-container").append("svg");
        //var diameter = 500;
        var div = d3.select("body").append("div")
            .attr("class", "tooltip")
            .style("opacity", 0);

        vis.attr("width", diameter)
            .attr("height", diameter);


        vis.selectAll("defs")
            .data(Object.keys(bt))
            .enter()
            .append("marker")
            .attr("id", function (d) {
            return ("lharrow_" + d);
        })
            .attr("viewBox", "-1 0 6 6")
            .attr("refX", -0.1)
            .attr("refY", 3)
            .attr("markerUnits", "strokeWidth")
            .attr("markerWidth", 1)
            .attr("markerHeight", 1)
            .attr("orient", "auto")
            .append("path")
            .attr("d", "M0,0 L-1,3 L0,6")
            .style("fill", function (d) {
            if (d == 'hilight') return 'yellow';
            return '#' + d;
            //return (bt_scale(d));
        });

        vis.selectAll("defs")
            .data(Object.keys(bt))
            .enter()
            .append("marker")
            .attr("id", function (d) {
            return ("rharrow_" + d);
        })
            .attr("viewBox", "0 0 6 6")
            .attr("refX", 0.1)
            .attr("refY", 3)
            .attr("markerUnits", "strokeWidth")
            .attr("markerWidth", 1)
            .attr("markerHeight", 1)
            .attr("orient", "auto")
            .append("path")
            .attr("d", "M0,0 L1,3 L0,6")
            .style("fill", function (d) {
            if (d == 'hilight') return 'yellow';
            return '#' + d;
            //return (bt_scale(d));
        });


        //add gene track

        vis.selectAll("svg")
            .data(genes)
            .enter()
            .append("path")
            .attr("d", function (d) {
            return (computeStrandPath(d.start, d.end, (diameter * 0.35) + (d.bumpLevel * 15), totalBP, "0,1"));
        })
            .attr("stroke-width", "10")
            .attr("transform", trans)
            .attr("fill", "none")
            .attr("stroke", function (d) {
            if (d.gene_name == $("#gene").val().toUpperCase()) {
                return ("yellow")
            } else {
                return '#' + d.color;
            }
        })
            .attr("marker-start", function (d) {
            var bt = d.color;
            if (d.gene_name == $("#gene").val().toUpperCase()) {
                bt = 'hilight';
            }

            if (d.strand == "-") return ("url(#lharrow_" + bt + ")");
        })
            .attr("marker-end", function (d) {
            var bt = d.color;
            if (d.gene_name == $("#gene").val().toUpperCase()) {
                bt = 'hilight';
            }
            if (d.strand == "+") return ("url(#rharrow_" + bt + ")");
        })
            .on("click", function (d) {
            console.log("Mouse Click !" + d.gene_name);
            $("#gene").val(d.gene_name);
            var gene = $("#gene").val().toUpperCase();
            div.transition()
                .duration(500)
                .style("opacity", 0);
            d3.selectAll("svg").remove();
            renderHic(gene, tissue, diameter, 1);
            return false;
        })
            .on("mouseover", function (d, i) {
            div.transition()
                .duration(200)
                .style("opacity", 0.9);
            div.html(d.gene_name + "</br>" + d.type + "</br>" + d.gene_id + "</br>" + numberWithCommas(parseInt(d.start) + parseInt(meta.ostart)) + "</br>" + numberWithCommas(parseInt(d.end) + parseInt(meta.ostart)))
                .style("left", (d3.event.pageX) + "px")
                .style("top", (d3.event.pageY - 28) + "px");

            d3.select(this).style("opacity", 0.3);
        })
            .on("mouseout", function (d) {
            div.transition()
                .duration(500)
                .style("opacity", 0);
            d3.select(this).style("opacity", 1);
        })


        // add hic links


        var color = d3.scale.linear()
            .domain([12, 20])
            .range(["blue", "red"]);

        vis.selectAll("svg")
        //.data(hics.filter(function(d){ return parseFloat(d.tissue[tissue]) >= score; }))
        .data(hics)
            .enter()
            .append("path")
            .attr("id", function (d, i) {
            return ('p' + i);
        })
            .attr("class", "interaction")
            .attr("d", function (d) {
            //return computePath(d.bSt,d.oeEnd,150,data.meta.rend - data.meta.rstart);
            return computePath(d.bSt + ((d.bEnd - d.bSt) / 2), d.oeSt + ((d.oeEnd - d.oeSt) / 2), diameter * 0.28, totalBP, diameter);
        })
            .attr("transform", trans)
            .attr("fill", "none")
            .attr("stroke", function (d) {
            return 'green';
            //return color(d.B23_CD4_Naive_D3_4)
            //return color(d.tissue[tissue]);
        })
            .attr("stroke-width", 4)
            .on("mouseover", function (d, i) {
            div.transition()
                .duration(200)
                .style("opacity", 0.9);
            //div.html("HicCup score " + parseFloat(d.tissue[tissue]).toFixed(2) + "</br>Bait: " + d.bSt + " " + d.bEnd + "</br>Target: " + d.oeSt + " " + d.oeEnd)
            var bStart = numberWithCommas(d.bSt + parseInt(meta.ostart));
            var bEnd = numberWithCommas(d.bEnd + parseInt(meta.ostart));
            var tStart = numberWithCommas(d.oeSt + parseInt(meta.ostart));
            var tEnd = numberWithCommas(d.oeEnd + parseInt(meta.ostart));
            div.html("Bait: " + bStart + '-' + bEnd + "</br>Target: " + tStart + '-' + tEnd)
                .style("left", (d3.event.pageX) + "px")
                .style("top", (d3.event.pageY - 28) + "px");

            // have a zscore issue
            // tried to use 'use' that doesn't work !
            // this approach is simple but does not help with 
            vis.select('#p' + i)
                .style("stroke", "yellow")
                .style("stroke-width", "10");

            // try sort approach from stackexchange http://stackoverflow.com/questions/13595175/updating-svg-element-z-index-with-d3

            vis.selectAll("path.interaction").sort(function (a, b) { // select the parent and sort the path's
                if (a.id != d.id) return -1; // a is not the hovered element, send "a" to the back
                else return 1; // a is the hovered element, bring "a" to the front
            });

            vis.append("path")
                .attr("class", "deleteMe")
                .attr("d", computeArcPath(d.oeSt, d.oeEnd, diameter * 0.28, diameter / 2, totalBP))
                .style("stroke-width", 1)
                .style("stroke", "red")
                .attr("transform", trans)
                .attr("fill", "none")

            vis.append("path")
                .attr("class", "deleteMe")
                .attr("d", computeArcPath(d.bSt, d.bEnd, diameter * 0.28, diameter / 2, totalBP))
                .style("stroke-width", 1)
                .style("stroke", "blue")
                .attr("transform", trans)
                .attr("fill", "none")


        })
            .on("mouseout", function (d, i) {
            div.transition()
                .duration(500)
                .style("opacity", 0);
            vis.select('#p' + i)
                .style("stroke", "green")
                .style("stroke-width", "4");
            vis.selectAll(".deleteMe").remove();

        })
            .on("click", function (d, i) {
            d3.select("#bait-svg").remove();
            d3.select("#target-svg").remove();
            $("#footer-bait").html('chr' + meta.rchr + ':' + numberWithCommas(d.bSt + parseInt(meta.ostart)) + '..' + numberWithCommas(d.bEnd + parseInt(meta.ostart)) + " (" + ((d.bEnd - d.bSt) / 1000).toFixed(2) + "KB)");
            $("#footer-target").html('chr' + meta.rchr + ':' + numberWithCommas(d.oeSt + parseInt(meta.ostart)) + '..' + numberWithCommas(d.oeEnd + parseInt(meta.ostart)) + " (" + ((d.oeEnd - d.oeSt) / 1000).toFixed(2) + "KB)");
            
            drawRegionPanel("bait", meta.rchr, (d.bSt + parseInt(meta.ostart)), (d.bEnd + parseInt(meta.ostart)), maxscore);
            drawRegionPanel("target", meta.rchr, (d.oeSt + parseInt(meta.ostart)), (d.oeEnd + parseInt(meta.ostart)), maxscore);
        });

        function log10(val) {
            return Math.log(val) / Math.LN10;
        }
        var thresh = -1 * log10(1e-1);
        //var thresh = -1 * log10(1e-5);
        var symb = d3.svg.symbol();
        symb.size(10);
        vis.selectAll("svg")
            .data(snps.filter(function (d) {
            return parseFloat(d.score) >= thresh;
        }))
            .enter()
            .append("path")
            .attr("transform", function (d) {
            //console.log(d);
            return (computePointPath(d.start, d.end, d.score, thresh, maxscore, diameter * 0.29, totalBP, diameter))
        })
            .attr("d", symb)
            .attr("stroke", function (d) {
            if (parseFloat(d.score) == maxscore) return "red";
            return "black";
        })
            .attr("fill", function (d) {
            if (parseFloat(d.score) == maxscore) return "red";
            return "black";
        })
            .on("mouseover", function (d, i) {
            div.transition()
                .duration(200)
                .style("opacity", 0.9);
            div.html(d.snp_name + "</br>" + d.score + "</br>" + numberWithCommas(parseInt(d.start) + parseInt(meta.rstart)) + "</br>")
                .style("left", (d3.event.pageX) + "px")
                .style("top", (d3.event.pageY - 28) + "px");

            d3.select(this).style("opacity", 0.3);
        })
            .on("mouseout", function (d) {
            div.transition()
                .duration(500)
                .style("opacity", 0);
            d3.select(this).style("opacity", 1);
        });

        vis.selectAll("svg")
            .data([1])
            .enter()

            .append("path")
        //.attr("d",computeArcPath(1,totalBP, diameter*0.2, diameter*0.1,  totalBP))
        .attr("d", computeStrandPath(1, totalBP, diameter * 0.35, totalBP, "0,0"))
            .style("stroke-width", 60)
            .style("stroke", "lightgrey")
            .attr("transform", trans)
            .attr("fill", "none");
        // end of JSON call   
    });
}

function drawRegionPanel(type, chr, start, end, maxscore) {	
	var region = chr+':'+start+'-'+end,
		data1 = [start, end],
		w = 700, h = 250, margin = 20, trackHeight = 100,
		formatxAxis = d3.format('.0f'),
		xRange = d3.scale.linear().domain([d3.min(data1), d3.max(data1)]).range([0 + margin, w - margin]),
		regionStart = d3.min(data1);
		
	//console.log(region);
		
	d3.json("/CHIC_DEMO2/cgi-bin/prototype.pl?region=" + region, function (error, data) {
			if (error) return console.warn(error);
			
			//console.log(data);
			adjustBump(data.genes, 0);	
			
			//Create the SVG Viewport
			var svgContainer = d3.select("#panel-"+type).append("svg:svg")
				.attr("width", w)
				.attr("height", h)
				.attr("id", type+"-svg");
	
			//Create the Axis
			var xAxis = d3.svg.axis()
				.scale(xRange)
				.ticks(4)
				.tickFormat(formatxAxis);
	
			//Create an SVG group Element for the Axis elements and call the xAxis function
			var xAxisGroup = svgContainer.append("svg:g")
				.attr("class", "x axis")
				.attr("transform", "translate(0," + (h - margin) + ")")
				.call(xAxis);	
				
			// TRACK 1 - SNPS
			var yRangeS = d3.scale.linear().domain([0, maxscore]).range([trackHeight,margin]);
			var yAxis = d3.svg.axis()
				.scale(yRangeS)
				.ticks(3)
				.tickFormat(d3.format('.0f'))
				.orient("left");
			
			var yAxisGroup = svgContainer.append("svg:g")
				.attr("class", "y axis")
				.attr("transform", "translate("+margin+",0)")
				.call(yAxis);
				
			var snp = svgContainer.append("g").attr("class", "track snps").selectAll(".snp")
				.data(data.snps)
				.enter().append("g")
				.attr("class", "snp");			
			
			snp.append("path")
				.attr("class", "marker")
				.attr("d", d3.svg.symbol().size(40))
				.attr("stroke", function (d) {
					if (parseFloat(d.score) == maxscore) return "red";
					if (parseFloat(d.score) >= 7.03) return "green";
					return "grey";
				})
				.attr("fill", function (d) {
					if (parseFloat(d.score) == maxscore) return "red";
					if (parseFloat(d.score) >= 7.03) return "green";
					return "grey";
				})
				.attr("transform", function (d) {
					return "translate(" + xRange(d.start + regionStart) + "," + yRangeS(d.score) + ")";
				});
			
			snp.append("text")
				.attr("y", function (d) {
				return yRangeS(d.score);
			})
				.attr("x", function (d) {
				return xRange(d.start + regionStart) - margin;
			})
				.attr("transform", function (d) {
				return "translate(" + (margin+5) + ",5)";
			})
				.text(function (d) {
				return d.snp_name;
			});
				
			// TRACK 2 - GENES
			var yRangeG = d3.scale.linear().domain([0, trackHeight]).range([margin, margin + trackHeight]);
			var geneTrackOffset = trackHeight + (2*margin);
			
			var line = d3.svg.line()
				.interpolate("linear")
				.x(function (d) { return xRange(d.x+regionStart); })
				.y(function (d) { return yRangeG(d.y+10); });
			
			var gene = svgContainer.append("g").attr("class", "track genes").selectAll(".gene")
				.data(data.genes)
				.enter().append("g")
				.attr("class", "gene");
				
			gene.append("path")
				.attr("class", "line")
				.attr("d", function (d) {
						return line([ { x: d.start, y: geneTrackOffset+(30 * d.bumpLevel)}, { x: d.end, y: geneTrackOffset+(30 * d.bumpLevel) }]);
				})
				.attr("stroke-width", "10")
				.attr("fill", "none")
				.attr("stroke", function (d) {
					return '#' + d.color;
				});
				
			gene.append("text")
				.attr("y", function(d) { return geneTrackOffset + yRangeG(30*d.bumpLevel); } )
				.attr("x", function (d) { return xRange(d.start + regionStart) - margin; })
				.attr("transform", function (d) {
					return "translate(" + margin + ",0)";
				})
				.text(function (d) {
					return d.gene_name;
				});
	});	
}

function resetPage(gene, tissue, breadcrumb) {
    d3.selectAll("svg").remove();
    $("#gene").val(gene);
    $("#page_header").html(gene + " in " + tissue.replace("_", " ") + " Tissues");
    $("#footer-bait").html("&nbsp;");
    $("#footer-target").html("&nbsp;");
    if (breadcrumb) $("#breadcrumb").append('<li id="BC-' + gene + '"><a href="#" onclick=\'javascript:d3.selectAll("svg").remove(); $("#BC-' + gene + '").remove(); renderHic("' + gene + '", $("input:radio[name=tissue]:checked").val(), 750, 1)\'>' + gene + '</a></li>');
}

$(document).ready(function () {
    $("#pushme").bind("click", function () {
        var tissue = $("input:radio[name=tissue]:checked").val();
        //var diameter = $("input:radio[name=di]:checked").val();
        var diameter = 750;
        var gene = $("#gene").val().toUpperCase();
        renderHic(gene, tissue, diameter, 1)
        return (false);
    });

    $("input:radio[name=tissue]").bind("click", function () {
        var tissue = $("input:radio[name=tissue]:checked").val();
        //var diameter = $("input:radio[name=di]:checked").val();
        var diameter = 750;
        var gene = $("#gene").val().toUpperCase();
        renderHic(gene, tissue, diameter, 0)
    });

    /*$("input:radio[name=di]").bind("click", function () {
        var tissue = $("input:radio[name=tissue]:checked").val();
        //var diameter = $("input:radio[name=di]:checked").val();
        var diameter = 750;
        var gene = $("#gene").val().toUpperCase();
        renderHic(gene, tissue, diameter, 0)
    });*/
});

var geneParam = getQueryVariable("gene");
if (geneParam == undefined){ geneParam = 'IL2RA'; }
renderHic(geneParam, 'CD4_Activated', 750, 1)