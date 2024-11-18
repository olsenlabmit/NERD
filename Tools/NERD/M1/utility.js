function matchAndReplace(patt, BigSmilesStr, func) {
    var n = 0;
    var match, indices = [];
    var repStrs = [];
    while (match = patt.exec(BigSmilesStr)) {
        n++;
        indices.push([match.index, match.index + match[0].length]);
        repStrs.push(func(n, match));
    }
    var i;
    for (i = n - 1; i >= 0; i--) {
        BigSmilesStr = BigSmilesStr.substr(0, indices[i][0]) + repStrs[i] + BigSmilesStr.substr(indices[i][1], BigSmilesStr.length + 1);
    }
    return BigSmilesStr;
}

function matchAndReplaceTwo(patt, BigSmilesStr, func) {
    var n = 0;
    var match, indices = [];
    var repStrs = [];
    var mass = [];
    var objects = $("#NERD").val() + "NEW\n" +
        "Temperature\n" + document.getElementById("inT").value + "\n" +
        "Number of KMC\n" + document.getElementById("inkmc").value + "\n" +
        "Conversion\n" + document.getElementById("inconv").value + "\n" +
        "Concentration\n" + document.getElementById("inconc").value + "\n" +
        "Number in Simulation\n" + document.getElementById("inmol").value + "\n"

    while (match = patt.exec(BigSmilesStr)) {
        n++;
        indices.push([match.index, match.index + match[0].length]);
        var list = func(n, match);
        repStrs.push(list[0]);
        mass.push([list[1], list[2], list[3], list[4], list[5]]);
    }
    var i;
    for (i = n - 1; i >= 0; i--) {
        BigSmilesStr = BigSmilesStr.substr(0, indices[i][0]) + repStrs[i] + BigSmilesStr.substr(indices[i][1], BigSmilesStr.length + 1);
    }
    objects += "Final BigSMILES\n" + BigSmilesStr + "\n" + 
                "Number of Objects\n" + n + "\n"
    for (i = 0; i < n; i++) {
        objects += "Object\n" + repStrs[i] + "\n" +
                    "Da\n" + mass[i][0] + "\n" +
                    "kDa\n" + mass[i][1] + "\n" +
                    "DP\n" + mass[i][2] + "\n" +
                    "Mo\n" + mass[i][3] + "\n" +
                    "b\n" + mass[i][4] + "\n";
    }
    objects += "\n"
    var list = [BigSmilesStr, objects];
    return list;
}

function escapeHTMLBrak(stringToShow, isNum) {
    var pattL = /&lt;/ig;
    var pattR = /&gt;/ig;
    stringToShow = stringToShow.replace(pattL, '<');
    stringToShow = stringToShow.replace(pattR, '>');
    return stringToShow
}


function correctBondDesc(BigSmilesStr) {
    BigSmilesStr = matchAndReplace(bondDescPatt, BigSmilesStr, function (x, match) {
        return match[1];
    });
    return BigSmilesStr;
}

function highlightStoObj(BigSmilesStr) {
    var patt = new RegExp(stoObjWidx.source, stoObjWidx.flags);
    BigSmilesStr = BigSmilesStr.replace(patt, '<span class="redSpan">$1</span>');
    patt = new RegExp(smilesCharPatt.source, smilesCharPatt.flags);
    BigSmilesStr = BigSmilesStr.replace(patt, '$2');
    patt = new RegExp(bondDescPatt.source, bondDescPatt.flags);
    BigSmilesStr = BigSmilesStr.replace(patt, '[$1]');
    // remove atom idx from smiles string
    patt = new RegExp(anyRegAtomWidx.source, anyRegAtomWidx.flags);
    BigSmilesStr = BigSmilesStr.replace(patt, '[$1]');
    return BigSmilesStr;
}


function getBigSmilesObjStr(localDivElement, isNum = false, useNumbered = true, get_bigsmiles = true) {
    //var localDivElement = element;
    var BigSmilesStr;
    var patt;
    if (isNum == false) {
        BigSmilesStr = localDivElement.children(".PolymerSMILES").html();
        patt = new RegExp(stoObjWidx.source, stoObjWidx.flags);
    } else {
        if (useNumbered == true) {
            BigSmilesStr = localDivElement.children(".PolymerNumberedSMILES").html();
        } else {
            BigSmilesStr = localDivElement.children(".PolymerSMILES").html();
        }
        var pattAtom = new RegExp(anyAtomWaddress.source, anyAtomWaddress.flags);
        var address = localDivElement.children(".title").children(".address").html();
        /*
            while(match = pattAtom.exec(BigSmilesStr)) {
            }
        */
        BigSmilesStr = BigSmilesStr.replace(pattAtom, "[$1:" + address + folder_div + "$2]");
        patt = new RegExp(stoObjWaddress.source, stoObjWaddress.flags);
    }

    // substitute the stochastic objects
    list = matchAndReplaceTwo(patt, BigSmilesStr, function (x, match) {
        return getStoObjStr(localDivElement.children(".StoObj").eq(x - 1), isNum, useNumbered);
    });
    BigSmilesStr = list[0];
    objects = list[1];
    //alert(BigSmilesStr);
    if (BigSmilesStr == '') {
        // complain if BigSmilesStr is empty
        var msg = "Warning: BigSMILES Unit    \"" + localDivElement.children(".title").children(".address").html() + "\" is empty!";
        $("#ErrMsg").val($("#ErrMsg").val() + msg + "\n");
    }
    if (get_bigsmiles)
        return BigSmilesStr;
    else
        return objects;
}

function getStoObjStr(localDivElement, isNum = false, useNumbered) {
    //alert('stoObjCalled');
    //var localDivElement = element;
    var RepUnitElement = localDivElement.children(".RepeatUnit");
    var EndGrpElement = localDivElement.children(".EndGroup");
    var BigSmilesStr = "{";
    // get left bonding descriptor
    var leftDesc = localDivElement.children(".DescText").children(".leftDesc").val();
    if (leftDesc == "") {
        //if(localDivElement.children(".DescText").children(".leftEnd")=='1'){
        // StoObj is left terminal, so possibly no left bond desc
        //} else {
        // StoObj is not left terminal, left bond desc must be present
        var msg = "Warning: Stochastic Object \"" + localDivElement.children(".title").children(".address").html() + "\" is missing bonding descriptor for " + escapeHTMLBrak(localDivElement.children(".DescText").children(".leftEndStr").html(), true);
        $("#ErrMsg").val($("#ErrMsg").val() + msg + "\n");
        //alert(msg);
        //}
        BigSmilesStr = BigSmilesStr + "[]"
    } else {
        BigSmilesStr = BigSmilesStr + "[" + leftDesc + "]"
    }
    // get all repeat units
    var nUnit = parseInt(RepUnitElement.children(".count").html());
    //alert(element.prop("tagName"));
    var i;
    for (i = 0; i < nUnit; i++) {
        if (i > 0) {
            BigSmilesStr = BigSmilesStr + ',';
        }
        BigSmilesStr = BigSmilesStr + getBigSmilesObjStr(RepUnitElement.children(".Polymer").eq(i), isNum, useNumbered);
    }

    if (nUnit == 0) {
        // complain if repeat unit list is empty
        var msg = "Warning: Stochastic Object \"" + localDivElement.children(".title").children(".address").html() + "\" has no repeat unit!";
        $("#ErrMsg").val($("#ErrMsg").val() + msg + "\n");

    }

    // get all end groups
    var nEnd = parseInt(EndGrpElement.children(".count").html());
    if (nEnd > 0) {
        BigSmilesStr = BigSmilesStr + ';';
    }
    for (i = 0; i < nEnd; i++) {
        if (i > 0) {
            BigSmilesStr = BigSmilesStr + ',';
        }
        BigSmilesStr = BigSmilesStr + getBigSmilesObjStr(EndGrpElement.children(".Polymer").eq(i), isNum, useNumbered);
    }
    // get right bonding descriptor
    var rightDesc = localDivElement.children(".DescText").children(".rightDesc").val();
    if (rightDesc == "") {
        if (parseInt(localDivElement.children(".DescText").children(".noEnds").html()) > 1) {
            var msg = "Warning: Stochastic Object \"" + localDivElement.children(".title").children(".address").html() + "\" is missing bonding descriptor for " + escapeHTMLBrak(localDivElement.children(".DescText").children(".rightEndStr").html(), true);
            $("#ErrMsg").val($("#ErrMsg").val() + msg + "\n");
        }

        //if(localDivElement.children(".DescText").children(".rightEnd")=='1'){
        // StoObj is right terminal, so possibly no right bond desc
        //} else {
        // StoObj is not right terminal, right bond desc must be present
        //var msg = "Warning: Stochastic Object \"" + localDivElement.children(".title").children(".address").html() + "\" is missing bonding descriptor!";
        //$("#ErrMsg").val( $("#ErrMsg").val() +  msg + "\n");
        //alert(msg);
        //}
        BigSmilesStr = BigSmilesStr + "[]"
    } else {
        BigSmilesStr = BigSmilesStr + "[" + rightDesc + "]"
    }
    var nTerminalDesc = parseInt(localDivElement.children(".DescText").children(".noEnds").html());
    if (nTerminalDesc > 2) {
        var msg = "Error: Degree of Stochastic Object \"" + localDivElement.children(".title").children(".address").html() + "\" larger than 2!";
        $("#ErrMsg").val($("#ErrMsg").val() + msg + "\n");
    }
    if (nEnd == 0 && nTerminalDesc < 2) {
        var msg = "Error: No implicit end group specified for Stochastic Object \"" + localDivElement.children(".title").children(".address").html() + "\", but only zero or one explicit end groups defined. ";
        $("#ErrMsg").val($("#ErrMsg").val() + msg + "\n");
    }

    BigSmilesStr = BigSmilesStr + '}';

    var inDa = localDivElement.children(".DescText").children(".inDa").val();
    var inkDa = localDivElement.children(".DescText").children(".inkDa").val();
    var inRU = localDivElement.children(".DescText").children(".inRU").val();  
    var inMo = localDivElement.children(".DescText").children(".inMo").val();
    var inb = localDivElement.children(".DescText").children(".inb").val();      
    return [BigSmilesStr, inDa, inkDa, inRU, inMo, inb];
}

function pad(n, width, z) {
    var padding = z || '0';
    n = n + '';
    return n.length >= width ? n : new Array(width - n.length + 1).join(padding) + n;
}

function numberAtoms(jmeStr, prefix) {
    var arrs = jmeStr.split("|");
    var atomCount = 0;
    var arr;
    for (var j = 0; j < arrs.length; j++) {
        arr = arrs[j].split(" ");
        var nAtom = parseInt(arr[0]);
        var nBond = parseInt(arr[1]);
        var idx;
        for (var i = 0; i < nAtom; i++) {
            atomCount++;
            idx = 2 + i * 3;
            //      arr[idx] = arr[idx] + ':'  + pad(atomCount, 2);

            var patt = new RegExp(bondDescPattinJME.source, bondDescPattinJME.flags);
            var match = patt.exec(arr[idx]);
            if (match == null) {
                //arr[idx] = arr[idx] + ':' + prefix + pad(atomCount, 2);
                arr[idx] = arr[idx] + ':' + pad(atomCount, 2);
            } else {
                arr[idx] = arr[idx];
            }

        }
        arrs[j] = arr.join(' ');
    }
    return arrs.join('|');
}

function clearPolymer(divElement) {
    divElement.children(".PolymerColoredSMILES").html("");
    divElement.children(".PolymerSMILES").html("");
    divElement.children(".PolymerNumberedSMILES").html("");
    divElement.children(".PolymerHidden").html("");
    divElement.children(".PolymerNumberedHidden").html("");
    divElement.children(".PolymerSVG").html("");
    divElement.children(".StoObj").remove();
    divElement.children(".Polymer").remove();
    //jsmeApplet.clear();
}

function getStoObjEnd(indexedStoObj, jmeStr) {
    var arrs = jmeStr.split("|");
    //  var atomCount=0;
    var arr;
    for (var j = 0; j < arrs.length; j++) {
        arr = arrs[j].split(" ");
        var nAtom = parseInt(arr[0]);
        var nBond = parseInt(arr[1]);
        var idx;
        for (var i = 0; i < nAtom; i++) {
            //    atomCount++;
            idx = 2 + i * 3;
            //arr[idx] = arr[idx] + ':' + prefix + pad(atomCount, 2);
            if (indexedStoObj == arr[idx]) {
                var returnArr = [];
                stoObjIdx = i + 1;
                for (var k = 0; k < nBond; k++) {
                    idx1 = 2 + nAtom * 3 + k * 3;
                    idx2 = 2 + nAtom * 3 + k * 3 + 1;
                    if (parseInt(arr[idx1]) == stoObjIdx) {
                        neighborIdx = parseInt(arr[idx2]);
                        returnArr.push(arr[(neighborIdx - 1) * 3 + 2]);
                    } else if (parseInt(arr[idx2]) == stoObjIdx) {
                        neighborIdx = parseInt(arr[idx1]);
                        returnArr.push(arr[(neighborIdx - 1) * 3 + 2]);
                    } else {
                        continue;
                    }
                }
                return returnArr;
            }
        }
        //    arrs[j] = arr.join(' ');
    }
    //return arrs.join('|');
}

function sortWithIndices(toSort) {
    for (var i = 0; i < toSort.length; i++) {
        toSort[i] = [toSort[i], i];
    }
    toSort.sort(function (left, right) {
        return left[0] < right[0] ? -1 : 1;
    });
    toSort.sortIndices = [];
    for (var j = 0; j < toSort.length; j++) {
        toSort.sortIndices.push(toSort[j][1]);
        toSort[j] = toSort[j][0];
    }
    return toSort;
}