#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3

// =============================================================================
// SMI_Layout.ipf - Layout Creation Functions
// =============================================================================
// Version 2.4.1 - 
// =============================================================================

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Static Constant kPtPerInch = 72
Static Constant kMmPerInch = 25.4
Static Constant kPageTitleHeight = 22
Static Constant kGraphLabelHeight = 18

// 
Static Constant kOutputGraph = 0    // Igor Graph ()
Static Constant kOutputPNG = 1      // PNG ()
Static Constant kOutputSVG = 2      // SVG ()

// 
// 
// 
// : Cell → Sample Average → Compare
Static StrConstant kGraphTypeOrder = "AvgAlignedTraj;AlignedTraj;Trajectory;MSD;StepHist;StepHeatmap;IntHist;LP;DstateDensity;ParticleDensity;MolDensDist;MolDensityHMM;DensityLogLog;XY_Density;DensityTime;AIC_Duration;Fraction_Duration;Tau_Duration;Duration;Onrate_Summary;Onrate;StateTransition;LowerBound;AverageHeatmap;AverageMSD;AverageStepHist;AverageIntHist;AverageLPHist;AverageStateFraction;AverageMolDensHist;AverageParticleVsMolDens;AverageOntime;AverageOnrate;AverageStateTrans;Compare_LowerBound;CompareMSD_D;CompareMSD_L;Compare_D;Compare_L;Compare_HMMP_All;Compare_HMMP;Compare_Tau_All;Compare_Tau;Compare_Fraction_All;Compare_Frac;Compare_OnRate_All;Compare_OnRate;Compare_Trans_All;Compare_Trans;Compare_Int_All;Compare_Int;Compare_LP_All;Compare_LP;Compare_MolDens;Compare_PDens;Compare_Area;Compare_NumPoints;Dstate"

//  → 
Static StrConstant kGraphTypeToButton = "AvgAlignedTraj:Avg Aligned;AlignedTraj:Aligned Traj;Trajectory:Trajectory;MSD:MSD;StepHist:Stepsize;StepHeatmap:Heatmap;IntHist:Intensity;LP:LP;DstateDensity:Mol Dens;ParticleDensity:Particle Dens;MolDensDist:Mol Dens Dist;MolDensityHMM:Mol Dens;DensityLogLog:Particle Dens;XY_Density:Particle Dens;DensityTime:Particle Dens;AIC_Duration:On-time;Fraction_Duration:On-time;Tau_Duration:On-time;Duration:On-time;Onrate_Summary:On-rate;Onrate:On-rate;StateTransition:State Trans;LowerBound:Lower Bound;AverageHeatmap:Avg Heatmap;AverageMSD:Avg MSD;AverageStepHist:Avg Step;AverageIntHist:Avg Int;AverageLPHist:Avg LP;AverageStateFraction:Avg State Frac;AverageMolDensHist:Avg MolDens;AverageParticleVsMolDens:Avg P vs M;AverageOntime:Avg Ontime;AverageOnrate:Avg Onrate;AverageStateTrans:Avg State Trans;Compare_LowerBound:Compare LB;CompareMSD_D:Compare MSD-D;CompareMSD_L:Compare MSD-L;Compare_D:Compare D;Compare_L:Compare L;Compare_HMMP_All:Compare Pop All;Compare_HMMP:Compare Pop;Compare_Tau_All:Compare Tau All;Compare_Tau:Compare Tau;Compare_Fraction_All:Compare Frac All;Compare_Frac:Compare Frac;Compare_OnRate_All:Compare OnRate All;Compare_OnRate:Compare OnRate;Compare_Trans_All:Compare Trans All;Compare_Trans:Compare Trans;Compare_Int_All:Compare Int All;Compare_Int:Compare Int;Compare_LP_All:Compare LP All;Compare_LP:Compare LP;Compare_MolDens:Compare MolDens;Compare_PDens:Compare PDens;Compare_Area:Compare Area;Compare_NumPoints:Compare NumPts;Dstate:D-state"

// =============================================================================
// : 
// =============================================================================
Function CreateAutoLayout(sampleName, pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	String sampleName
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	// CreateAutoLayoutWithList
	String graphList = GetAllGraphWindows()
	return CreateAutoLayoutWithList(sampleName, pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

// =============================================================================
// graphList
// =============================================================================
Function CreateAutoLayoutWithList(sampleName, pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
	String sampleName
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	String graphList
	
	// No Image title
	NVAR/Z noTitle = root:LayoutNoTitle
	Variable hideLabel = 0
	if(NVAR_Exists(noTitle))
		hideLabel = noTitle
	endif
	
	// MolDensImg_*
	String firstGraphInList = StringFromList(0, graphList)
	Variable isMolDensImg = StringMatch(firstGraphInList, "MolDensImg_*")
	
	// 
	Variable offset_pt = offset_mm / kMmPerInch * kPtPerInch
	Variable gap_pt = gap_mm / kMmPerInch * kPtPerInch
	Variable pageW_pt = pageW_inch * kPtPerInch
	Variable pageH_pt = pageH_inch * kPtPerInch
	
	// 
	Variable numGraphs = ItemsInList(graphList)
	
	if(numGraphs == 0)
		DoAlert 0, "No graph windows found."
		return 0
	endif
	
	Printf "Found %d graphs (output mode: %d)\r", numGraphs, outputMode
	
	// 
	Variable contentW = pageW_pt - 2 * offset_pt
	Variable contentH = pageH_pt - 2 * offset_pt - kPageTitleHeight
	
	// 
	Variable cellW = (contentW - (divW - 1) * gap_pt) / divW
	Variable cellH = (contentH - (divH - 1) * gap_pt) / divH
	
	// MolDensImg_*noTitle
	Variable maxGraphW = cellW
	Variable maxGraphH
	Variable effectiveLabelHeight
	if(isMolDensImg && hideLabel)
		maxGraphH = cellH
		effectiveLabelHeight = 0
	else
		maxGraphH = cellH - kGraphLabelHeight
		effectiveLabelHeight = kGraphLabelHeight
	endif
	
	// 
	Variable graphsPerPage = divW * divH
	
	// 
	String typeGroups = GroupGraphsByType(graphList)
	Variable numGroups = ItemsInList(typeGroups, "|")
	
	// 
	String layoutName = GenerateUniqueLayoutName(sampleName)
	
	// PNG/SVG
	if(outputMode != kOutputGraph)
		NewPath/O/Q SMI_TempPath, SpecialDirPath("Temporary", 0, 0, 0)
	endif
	
	// MolDensImg_*
	String cmd, layoutWindowTitle
	if(isMolDensImg)
		layoutWindowTitle = "Molecular Density Image"
	else
		layoutWindowTitle = "SMI Analysis Layout"
	endif
	sprintf cmd, "NewLayout/K=1/W=(50,50,%d,%d)/N=%s as \"%s\"", 50+round(pageW_pt/2), 50+round(pageH_pt/2), layoutName, layoutWindowTitle
	Execute cmd
	
	sprintf cmd, "LayoutPageAction/W=%s size(-1)=(%d, %d), margins(-1)=(18, 18, 18, 18)", layoutName, round(pageW_pt), round(pageH_pt)
	Execute cmd
	
	// 
	Variable globalGraphIdx = 0
	Variable pageNum = 0
	Variable groupIdx, graphIdxInGroup, row, col
	Variable numInGroup, posOnPage, k
	Variable cellLeft, cellTop, graphW, graphH, posX, posY
	String graphName, graphLabel, currentGroup, groupType, pageGraphList, pageTitle
	
	for(groupIdx = 0; groupIdx < numGroups; groupIdx += 1)
		currentGroup = StringFromList(groupIdx, typeGroups, "|")
		numInGroup = ItemsInList(currentGroup)
		
		// 
		if(numInGroup > 0)
			groupType = GetButtonNameFromGraph(StringFromList(0, currentGroup))
		else
			groupType = "Unknown"
		endif
		
		// 
		graphIdxInGroup = 0
		do
			// 
			posOnPage = mod(graphIdxInGroup, graphsPerPage)
			
			if(posOnPage == 0)
				// 
				if(pageNum > 0)
					sprintf cmd, "LayoutPageAction/W=%s appendPage", layoutName
					Execute cmd
				endif
				
				// 
				pageGraphList = ""
				for(k = 0; k < min(graphsPerPage, numInGroup - graphIdxInGroup); k += 1)
					pageGraphList += StringFromList(graphIdxInGroup + k, currentGroup) + ";"
				endfor
				
				// 
				pageTitle = GeneratePageTitle(pageGraphList, pageNum+1, -1)
				AddPageTitle(layoutName, pageTitle, pageW_pt, offset_pt, pageNum+1)
				
				pageNum += 1
			endif
			
			// 
			graphName = StringFromList(graphIdxInGroup, currentGroup)
			
			col = mod(posOnPage, divW)
			row = floor(posOnPage / divW)
			
			// 
			cellLeft = offset_pt + col * (cellW + gap_pt)
			cellTop = offset_pt + kPageTitleHeight + row * (cellH + gap_pt)
			
			// 
			GetFittedGraphSize(graphName, maxGraphW, maxGraphH, graphW, graphH)
			
			// 
			posX = cellLeft + cellW - graphW
			posY = cellTop
			
			// SampleName
			graphLabel = GetGraphEmbeddedTitle(graphName)
			
			// 
			if(outputMode == kOutputGraph)
				AppendAndResizeGraph(layoutName, graphName, posX, posY + effectiveLabelHeight, graphW, graphH)
			elseif(outputMode == kOutputPNG)
				AppendGraphAsPNG(layoutName, graphName, posX, posY + effectiveLabelHeight, graphW, graphH, globalGraphIdx)
			elseif(outputMode == kOutputSVG)
				AppendGraphAsSVG(layoutName, graphName, posX, posY + effectiveLabelHeight, graphW, graphH, globalGraphIdx)
			endif
			
			// MolDensImg_*hideLabel
			if(!(isMolDensImg && hideLabel))
				AddGraphLabel(layoutName, graphLabel, posX, posY, graphW, pageNum, globalGraphIdx)
			endif
			
			graphIdxInGroup += 1
			globalGraphIdx += 1
		while(graphIdxInGroup < numInGroup)
	endfor
	
	// 
	sprintf cmd, "ModifyLayout/W=%s mag=0.75", layoutName
	Execute cmd
	
	// 
	if(outputMode != kOutputGraph)
		KillPath/Z SMI_TempPath
	endif
	
	Printf "Layout '%s' created with %d pages\r", layoutName, pageNum
	return pageNum
End

// -----------------------------------------------------------------------------
// FolderName
// -----------------------------------------------------------------------------
Function/S GroupGraphsByType(graphList)
	String graphList
	
	String typeOrder = kGraphTypeOrder
	Variable numTypes = ItemsInList(typeOrder)
	Variable numGraphs = ItemsInList(graphList)
	
	String result = ""
	String usedGraphs = ""
	
	Variable i, j
	String typeName, graphName, currentGroup
	
	// 
	for(i = 0; i < numTypes; i += 1)
		typeName = StringFromList(i, typeOrder)
		currentGroup = ""
		
		// 
		for(j = 0; j < numGraphs; j += 1)
			graphName = StringFromList(j, graphList)
			if(WhichListItem(graphName, usedGraphs) >= 0)
				continue  // 
			endif
			
			if(StringMatch(graphName, "*" + typeName + "*"))
				currentGroup += graphName + ";"
				usedGraphs += graphName + ";"
			endif
		endfor
		
		// FolderName
		if(strlen(currentGroup) > 0)
			currentGroup = SortGraphsByFolderName(currentGroup)
			
			if(strlen(result) > 0)
				result += "|"
			endif
			result += currentGroup
		endif
	endfor
	
	// 
	currentGroup = ""
	for(j = 0; j < numGraphs; j += 1)
		graphName = StringFromList(j, graphList)
		if(WhichListItem(graphName, usedGraphs) < 0)
			currentGroup += graphName + ";"
		endif
	endfor
	
	if(strlen(currentGroup) > 0)
		currentGroup = SortGraphsByFolderName(currentGroup)
		if(strlen(result) > 0)
			result += "|"
		endif
		result += currentGroup
	endif
	
	return result
End

// -----------------------------------------------------------------------------
// FolderName
// Compare/Average(S0,S1)(C1,C2)
// _All 
// -----------------------------------------------------------------------------
Function/S SortGraphsByFolderName(graphList)
	String graphList
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs <= 1)
		return graphList
	endif
	
	// 
	Make/FREE/T/N=(numGraphs) graphNames
	Make/FREE/N=(numGraphs) sortKeys
	
	Variable i, numVal
	String graphName, folderName, sortedList
	
	for(i = 0; i < numGraphs; i += 1)
		graphName = StringFromList(i, graphList)
		graphNames[i] = graphName
		
		// _All 
		if(StringMatch(graphName, "*_All") || StringMatch(graphName, "*_All_*"))
			numVal = 9999
		// Compare/Average
		elseif(StringMatch(graphName, "Compare_*") || StringMatch(graphName, "Average*") || StringMatch(graphName, "CompareMSD*"))
			numVal = ExtractStateOrComponentNumber(graphName)
			if(numVal < 0)
				numVal = 5000 + i  // 
			endif
		else
			// FolderName
			folderName = GetFolderNameFromGraph(graphName)
			numVal = ExtractNumberFromString(folderName)
			
			// 
			if(numVal < 0)
				numVal = ExtractNumberFromString(graphName)
			endif
			
			// 
			if(numVal < 0)
				numVal = 10000 + i
			endif
		endif
		
		sortKeys[i] = numVal
	endfor
	
	// 
	Sort sortKeys, graphNames, sortKeys
	
	// 
	sortedList = ""
	for(i = 0; i < numGraphs; i += 1)
		sortedList += graphNames[i] + ";"
	endfor
	
	return sortedList
End

// -----------------------------------------------------------------------------
// Compare/Average(S0,S1)(C1,C2)
// -----------------------------------------------------------------------------
Function ExtractStateOrComponentNumber(graphName)
	String graphName
	
	Variable len = strlen(graphName)
	Variable i, numStart, numEnd, j, charCode
	String ch, numStr
	
	// _S0, _S1, _C1, _C2 
	//  _S  _C 
	for(i = len - 1; i >= 1; i -= 1)
		ch = graphName[i-1, i]
		if(StringMatch(ch, "_S") || StringMatch(ch, "_C"))
			// 
			numStart = i + 1
			if(numStart < len)
				numStr = ""
				for(j = numStart; j < len; j += 1)
					charCode = char2num(graphName[j, j])
					if(charCode >= 48 && charCode <= 57)  // '0'-'9'
						numStr += graphName[j, j]
					else
						break
					endif
				endfor
				if(strlen(numStr) > 0)
					return str2num(numStr)
				endif
			endif
		endif
	endfor
	
	// _dt1_S0 
	return ExtractNumberFromString(graphName)
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function ExtractNumberFromString(str)
	String str
	
	Variable len = strlen(str)
	Variable i, charCode, foundNumber
	String numStr, ch
	
	numStr = ""
	foundNumber = 0
	
	// 
	for(i = len - 1; i >= 0; i -= 1)
		ch = str[i, i]
		charCode = char2num(ch)
		
		if(charCode >= 48 && charCode <= 57)  // '0'-'9'
			numStr = ch + numStr
			foundNumber = 1
		elseif(foundNumber)
			// 
			break
		endif
	endfor
	
	if(strlen(numStr) > 0)
		return str2num(numStr)
	endif
	
	return -1
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function/S GenerateUniqueLayoutName(sampleName)
	String sampleName
	
	String baseName = "SMI_Layout"
	if(strlen(sampleName) > 0 && !StringMatch(sampleName, "(No samples loaded)"))
		baseName = "Layout_" + CleanupName(sampleName, 0)
	endif
	
	// 
	String layoutList = WinList("*", ";", "WIN:4")
	Variable counter = 0
	String testName = baseName
	
	do
		if(WhichListItem(testName, layoutList) < 0)
			// 
			return testName
		endif
		counter += 1
		sprintf testName, "%s_%d", baseName, counter
	while(counter < 100)
	
	return testName
End

// -----------------------------------------------------------------------------
//  + 
// numPages = -1 
// -----------------------------------------------------------------------------
Function/S GeneratePageTitle(pageGraphList, pageNum, numPages)
	String pageGraphList
	Variable pageNum, numPages
	
	Variable numGraphs = ItemsInList(pageGraphList)
	if(numGraphs == 0)
		return "SMI Analysis"
	endif
	
	// MolDensImg_*
	String firstGraph = StringFromList(0, pageGraphList)
	if(StringMatch(firstGraph, "MolDensImg_*"))
		String molDensTitle
		sprintf molDensTitle, "Molecular Density Image - Page %d", pageNum
		return molDensTitle
	endif
	
	// : 
	Variable i
	String graphName, oldLabel
	
	// 
	String commonFolder = ""
	String firstLabel = GetGraphLabel(StringFromList(0, pageGraphList))
	Variable colonPos = strsearch(firstLabel, ": ", 0)
	if(colonPos > 0)
		commonFolder = firstLabel[colonPos+2, strlen(firstLabel)-1]
		
		// 
		for(i = 1; i < numGraphs; i += 1)
			graphName = StringFromList(i, pageGraphList)
			String thisLabel = GetGraphLabel(graphName)
			Variable thisColonPos = strsearch(thisLabel, ": ", 0)
			if(thisColonPos > 0)
				String thisFolder = thisLabel[thisColonPos+2, strlen(thisLabel)-1]
				if(!StringMatch(thisFolder, commonFolder))
					commonFolder = ""  // 
					break
				endif
			endif
		endfor
	endif
	
	// 
	String typeList = ""
	for(i = 0; i < numGraphs; i += 1)
		graphName = StringFromList(i, pageGraphList)
		String thisLbl = GetGraphLabel(graphName)
		colonPos = strsearch(thisLbl, ": ", 0)
		String typePart
		if(colonPos > 0)
			typePart = thisLbl[0, colonPos-1]
		else
			typePart = thisLbl
		endif
		// 
		if(WhichListItem(typePart, typeList, ", ") < 0)
			if(strlen(typeList) > 0)
				typeList += ", "
			endif
			typeList += typePart
		endif
	endfor
	
	// 
	String pageTitle
	if(strlen(commonFolder) > 0)
		pageTitle = commonFolder + " (" + typeList + ")"
	else
		pageTitle = typeList
	endif
	
	// numPages > 0 
	if(numPages > 1)
		sprintf pageTitle, "%s - Page %d/%d", pageTitle, pageNum, numPages
	endif
	
	return pageTitle
End

// -----------------------------------------------------------------------------
// waveFolderName
// -----------------------------------------------------------------------------
Function/S GetGraphEmbeddedTitle(graphName)
	String graphName
	
	// MolDensImg_*ImageSampleName + Sn
	if(StringMatch(graphName, "MolDensImg_*"))
		String imageName = GetMolDensImageName(graphName)
		String molDensTitle = ExtractMolDensSampleAndState(imageName)
		if(strlen(molDensTitle) > 0)
			return molDensTitle
		endif
	endif
	
	// Compare/Average
	if(StringMatch(graphName, "Compare_*") || StringMatch(graphName, "Average*"))
		String windowTitle = GetWindowTitleForLabel(graphName)
		if(strlen(windowTitle) > 0)
			return windowTitle
		endif
	endif
	
	// FolderName
	String folderName = GetFolderNameFromGraph(graphName)
	
	if(strlen(folderName) > 0)
		return folderName
	endif
	
	// 
	return GetButtonNameFromGraph(graphName)
End

// -----------------------------------------------------------------------------
// MolDensImg_*Image
// -----------------------------------------------------------------------------
Function/S GetMolDensImageName(graphName)
	String graphName
	
	String imageList = ImageNameList(graphName, ";")
	if(strlen(imageList) > 0)
		return StringFromList(0, imageList)
	endif
	return ""
End

// -----------------------------------------------------------------------------
// MolDensImg ImageSampleName + Sn
// Image_S{state}_C{cond}_{sampleName} → {sampleName} S{state}
// -----------------------------------------------------------------------------
Function/S ExtractMolDensSampleAndState(imageName)
	String imageName
	
	if(strlen(imageName) == 0)
		return ""
	endif
	
	// "Image_" 
	String remaining = ReplaceString("Image_", imageName, "")
	
	// S{state} 
	Variable sPos = strsearch(remaining, "S", 0)
	Variable cPos = strsearch(remaining, "_C", 0)
	if(sPos < 0 || cPos < 0)
		return imageName
	endif
	
	String stateStr = remaining[sPos, cPos-1]
	
	// C{cond}_ SampleName
	Variable underscorePos = strsearch(remaining, "_", cPos+2)
	if(underscorePos < 0)
		return imageName
	endif
	
	String sampleName = remaining[underscorePos+1, strlen(remaining)-1]
	
	return sampleName + " " + stateStr
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function/S GetWindowTitleForLabel(graphName)
	String graphName
	
	DoWindow $graphName
	if(V_Flag == 0)
		return ""
	endif
	
	// 
	GetWindow $graphName, title
	String windowTitle = S_value
	
	// 
	if(strlen(windowTitle) > 40)
		// 
		Variable colonPos = strsearch(windowTitle, ":", 0)
		Variable parenPos = strsearch(windowTitle, "(", 0)
		
		if(parenPos > 0 && parenPos < 40)
			windowTitle = windowTitle[0, parenPos-1]
		elseif(colonPos > 0 && colonPos < 40)
			windowTitle = windowTitle[0, colonPos-1]
		else
			windowTitle = windowTitle[0, 37] + "..."
		endif
	endif
	
	// 
	windowTitle = SMI_TrimString(windowTitle)
	
	return windowTitle
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function/S SMI_TrimString(str)
	String str
	
	Variable len = strlen(str)
	Variable startPos = 0
	Variable endPos = len - 1
	
	// 
	do
		if(startPos >= len)
			return ""
		endif
		if(char2num(str[startPos, startPos]) != 32)  // 
			break
		endif
		startPos += 1
	while(1)
	
	// 
	do
		if(endPos < startPos)
			return ""
		endif
		if(char2num(str[endPos, endPos]) != 32)  // 
			break
		endif
		endPos -= 1
	while(1)
	
	return str[startPos, endPos]
End

// -----------------------------------------------------------------------------
// waveFolderName
// wave: root:SampleName:FolderName:wavename
// -----------------------------------------------------------------------------
Function/S GetFolderNameFromGraph(graphName)
	String graphName
	
	DoWindow $graphName
	if(V_Flag == 0)
		return ""
	endif
	
	// 
	String traceList = TraceNameList(graphName, ";", 1)
	Variable numTraces = ItemsInList(traceList)
	
	if(numTraces == 0)
		return ""
	endif
	
	// wave
	String firstTrace = StringFromList(0, traceList)
	Wave/Z w = TraceNameToWaveRef(graphName, firstTrace)
	
	if(!WaveExists(w))
		return ""
	endif
	
	// wavewave
	String folderPath = GetWavesDataFolder(w, 1)
	
	// FolderName (root:SampleName:FolderName: )
	// root:
	if(StringMatch(folderPath, "root:*"))
		folderPath = folderPath[5, strlen(folderPath)-1]
	endif
	
	// 
	if(StringMatch(folderPath[strlen(folderPath)-1, strlen(folderPath)-1], ":"))
		folderPath = folderPath[0, strlen(folderPath)-2]
	endif
	
	// 2FolderName
	// SampleName:FolderName  FolderName 
	Variable colonPos = strsearch(folderPath, ":", 0)
	if(colonPos > 0)
		String folderName = folderPath[colonPos+1, strlen(folderPath)-1]
		
		// root:Sample:Folder:SubFolder
		Variable nextColon = strsearch(folderName, ":", 0)
		if(nextColon > 0)
			folderName = folderName[0, nextColon-1]
		endif
		
		// 
		if(strlen(folderName) > 30)
			folderName = folderName[0, 27] + "..."
		endif
		
		return folderName
	endif
	
	// SampleName
	return folderPath
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function/S GetButtonNameFromGraph(graphName)
	String graphName
	
	String typeOrder = kGraphTypeOrder
	String typeToButton = kGraphTypeToButton
	Variable numTypes = ItemsInList(typeOrder)
	Variable i
	String typeName, buttonName
	
	for(i = 0; i < numTypes; i += 1)
		typeName = StringFromList(i, typeOrder)
		if(StringMatch(graphName, "*" + typeName + "*"))
			// 
			buttonName = StringByKey(typeName, typeToButton)
			if(strlen(buttonName) > 0)
				return buttonName
			endif
		endif
	endfor
	
	// 
	if(strlen(graphName) > 25)
		return graphName[0,22] + "..."
	endif
	return graphName
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function GetFittedGraphSize(graphName, maxW, maxH, outW, outH)
	String graphName
	Variable maxW, maxH
	Variable &outW, &outH
	
	// 
	GetWindow $graphName, wsize
	Variable winW = V_right - V_left
	Variable winH = V_bottom - V_top
	
	if(winW <= 0 || winH <= 0)
		outW = maxW
		outH = maxH
		return 0
	endif
	
	// 
	Variable aspect = winW / winH
	
	// 
	Variable scaleW = maxW / winW
	Variable scaleH = maxH / winH
	Variable scaleFactor = min(scaleW, scaleH)
	
	outW = winW * scaleFactor
	outH = winH * scaleFactor
	
	return 1
End

// -----------------------------------------------------------------------------
// PNG
// -----------------------------------------------------------------------------
Function AppendGraphAsPNG(layoutName, graphName, left, top, targetW, targetH, idx)
	String layoutName, graphName
	Variable left, top, targetW, targetH, idx
	
	DoWindow $graphName
	if(V_Flag == 0)
		return 0
	endif
	
	String cmd
	String pictName = "SMI_PNG_" + num2str(idx)
	String tempFileName = pictName + ".png"
	
	// PNG288 DPI
	sprintf cmd, "SavePICT/O/E=-5/B=288/WIN=%s/P=SMI_TempPath as \"%s\"", graphName, tempFileName
	Execute cmd
	
	// 
	sprintf cmd, "LoadPICT/O/Q/P=SMI_TempPath \"%s\", %s", tempFileName, pictName
	Execute cmd
	
	// 
	sprintf cmd, "AppendLayoutObject/W=%s/F=0/T=0 picture %s", layoutName, pictName
	Execute cmd
	
	// 
	sprintf cmd, "ModifyLayout/W=%s left(%s)=%d, top(%s)=%d, width(%s)=%d, height(%s)=%d", layoutName, pictName, round(left), pictName, round(top), pictName, round(targetW), pictName, round(targetH)
	Execute cmd
	
	return 1
End

// -----------------------------------------------------------------------------
// SVG
// -----------------------------------------------------------------------------
Function AppendGraphAsSVG(layoutName, graphName, left, top, targetW, targetH, idx)
	String layoutName, graphName
	Variable left, top, targetW, targetH, idx
	
	DoWindow $graphName
	if(V_Flag == 0)
		return 0
	endif
	
	String cmd
	String pictName = "SMI_SVG_" + num2str(idx)
	String tempFileName = pictName + ".svg"
	
	// SVG
	sprintf cmd, "SavePICT/O/E=-9/WIN=%s/P=SMI_TempPath as \"%s\"", graphName, tempFileName
	Execute cmd
	
	// 
	sprintf cmd, "LoadPICT/O/Q/P=SMI_TempPath \"%s\", %s", tempFileName, pictName
	Execute cmd
	
	// 
	sprintf cmd, "AppendLayoutObject/W=%s/F=0/T=0 picture %s", layoutName, pictName
	Execute cmd
	
	// 
	sprintf cmd, "ModifyLayout/W=%s left(%s)=%d, top(%s)=%d, width(%s)=%d, height(%s)=%d", layoutName, pictName, round(left), pictName, round(top), pictName, round(targetW), pictName, round(targetH)
	Execute cmd
	
	return 1
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function AppendAndResizeGraph(layoutName, graphName, left, top, targetW, targetH)
	String layoutName, graphName
	Variable left, top, targetW, targetH
	
	DoWindow $graphName
	if(V_Flag == 0)
		return 0
	endif
	
	String cmd
	
	sprintf cmd, "AppendLayoutObject/W=%s/F=0/T=0 graph %s", layoutName, graphName
	Execute cmd
	
	sprintf cmd, "ModifyLayout/W=%s left(%s)=%d, top(%s)=%d, width(%s)=%d, height(%s)=%d", layoutName, graphName, round(left), graphName, round(top), graphName, round(targetW), graphName, round(targetH)
	Execute cmd
	
	return 1
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function/S GetAllGraphWindows()
	String allGraphs = WinList("*", ";", "WIN:1")
	Variable numGraphs = ItemsInList(allGraphs)
	Variable i
	String graphName
	
	// 
	String filteredList = ""
	for(i = 0; i < numGraphs; i += 1)
		graphName = StringFromList(i, allGraphs)
		
		if(StringMatch(graphName, "SMI_*Panel*"))
			continue
		endif
		if(StringMatch(graphName, "Panel*"))
			continue
		endif
		
		filteredList += graphName + ";"
	endfor
	
	// 
	String sortedList = SortGraphsByType(filteredList)
	
	return sortedList
End

// -----------------------------------------------------------------------------
// FolderName
// -----------------------------------------------------------------------------
Function/S SortGraphsByType(graphList)
	String graphList
	
	String typeOrder = kGraphTypeOrder
	Variable numTypes = ItemsInList(typeOrder)
	Variable numGraphs = ItemsInList(graphList)
	
	String sortedList = ""
	String usedGraphs = ""
	
	Variable i, j
	String typeName, graphName, typeGroup
	
	// 
	for(i = 0; i < numTypes; i += 1)
		typeName = StringFromList(i, typeOrder)
		typeGroup = ""
		
		// 
		for(j = 0; j < numGraphs; j += 1)
			graphName = StringFromList(j, graphList)
			if(WhichListItem(graphName, usedGraphs) >= 0)
				continue
			endif
			
			if(StringMatch(graphName, "*" + typeName + "*"))
				typeGroup += graphName + ";"
				usedGraphs += graphName + ";"
			endif
		endfor
		
		// FolderName
		if(strlen(typeGroup) > 0)
			typeGroup = SortGraphsByFolderName(typeGroup)
			sortedList += typeGroup
		endif
	endfor
	
	// 
	typeGroup = ""
	for(j = 0; j < numGraphs; j += 1)
		graphName = StringFromList(j, graphList)
		if(WhichListItem(graphName, usedGraphs) < 0)
			typeGroup += graphName + ";"
		endif
	endfor
	
	if(strlen(typeGroup) > 0)
		typeGroup = SortGraphsByFolderName(typeGroup)
		sortedList += typeGroup
	endif
	
	return sortedList
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function/S GetGraphLabel(graphName)
	String graphName
	
	String typeOrder = kGraphTypeOrder
	String typeToButton = kGraphTypeToButton
	Variable numTypes = ItemsInList(typeOrder)
	Variable i
	String typeName, buttonName
	
	for(i = 0; i < numTypes; i += 1)
		typeName = StringFromList(i, typeOrder)
		if(StringMatch(graphName, "*" + typeName + "*"))
			buttonName = StringByKey(typeName, typeToButton)
			if(strlen(buttonName) > 0)
				String folderPart = ExtractFolderFromGraphName(graphName, typeName)
				if(strlen(folderPart) > 0)
					return buttonName + ": " + folderPart
				else
					return buttonName
				endif
			endif
		endif
	endfor
	
	if(strlen(graphName) > 25)
		return graphName[0,22] + "..."
	endif
	return graphName
End

// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
Function/S ExtractFolderFromGraphName(graphName, typeName)
	String graphName, typeName
	
	Variable typePos = strsearch(graphName, typeName, 0, 2)
	String folderPart = ""
	
	if(typePos >= 0)
		if(typePos > 0)
			folderPart = graphName[0, typePos-2]
		else
			Variable afterType = typePos + strlen(typeName)
			if(afterType < strlen(graphName) - 1)
				folderPart = graphName[afterType+1, strlen(graphName)-1]
			endif
		endif
	endif
	
	folderPart = ReplaceString("_", folderPart, " ")
	
	if(strlen(folderPart) > 18)
		folderPart = folderPart[0,15] + "..."
	endif
	
	return folderPart
End

// -----------------------------------------------------------------------------
// Arial 18pt, Bold
// -----------------------------------------------------------------------------
Function AddPageTitle(layoutName, titleText, pageW, offset, pageNum)
	String layoutName, titleText
	Variable pageW, offset, pageNum
	
	String textName = "PageTitle_" + num2str(pageNum)
	String cmd
	
	sprintf cmd, "TextBox/W=%s/C/N=%s/F=0/A=LT/X=0/Y=0/E=2 \"\\F'Arial'\\Z18\\f01%s\"", layoutName, textName, titleText
	Execute cmd
	
	sprintf cmd, "ModifyLayout/W=%s frame(%s)=0, left(%s)=%g, top(%s)=%g", layoutName, textName, textName, offset, textName, offset
	Execute cmd
End

// -----------------------------------------------------------------------------
// Arial 14pt
// -----------------------------------------------------------------------------
Function AddGraphLabel(layoutName, labelText, left, top, graphWidth, pageNum, idx)
	String layoutName, labelText
	Variable left, top, graphWidth, pageNum, idx
	
	String textName = "Label_" + num2str(pageNum) + "_" + num2str(idx)
	String cmd
	
	sprintf cmd, "TextBox/W=%s/C/N=%s/F=0/A=LT/X=0/Y=0/E=2 \"\\F'Arial'\\Z14%s\"", layoutName, textName, labelText
	Execute cmd
	
	sprintf cmd, "ModifyLayout/W=%s frame(%s)=0, left(%s)=%g, top(%s)=%g", layoutName, textName, textName, left, textName, top
	Execute cmd
End

// =============================================================================
// 
// =============================================================================

Function MmToPoints(mm)
	Variable mm
	return mm / kMmPerInch * kPtPerInch
End

Function InchToPoints(inch)
	Variable inch
	return inch * kPtPerInch
End

// =============================================================================
// Colocalization Layout Functions
// =============================================================================

// Colocalization
// graphType: "Batch", "Average", "Compare", "All"
Function/S GetColocalizationGraphList(graphType)
	String graphType
	
	String graphList = ""
	String allWindows = WinList("*", ";", "WIN:1")  // All graphs
	Variable numWindows = ItemsInList(allWindows)
	Variable i
	String winName
	
	for(i = 0; i < numWindows; i += 1)
		winName = StringFromList(i, allWindows)
		
		// 
		Variable include = 0
		
		if(StringMatch(graphType, "Batch") || StringMatch(graphType, "All"))
			// Batch Analysis graphs
			if(StringMatch(winName, "*_ColTraj*"))  // Colocalization Trajectory
				include = 1
			endif
			if(StringMatch(winName, "Dist_*_C*E*"))  // Distance histogram
				include = 1
			endif
			if(StringMatch(winName, "IntHist_*_C*E*"))  // Intensity histogram
				include = 1
			endif
			if(StringMatch(winName, "StepHist_*_C*E*"))  // Stepsize histogram
				include = 1
			endif
			if(StringMatch(winName, "MSD_*_C*E*"))  // MSD
				include = 1
			endif
			if(StringMatch(winName, "Duration_*_C*E*"))  // On-time
				include = 1
			endif
			if(StringMatch(winName, "Onrate_*_C*E*"))  // On-rate
				include = 1
			endif
		endif
		
		if(StringMatch(graphType, "Average") || StringMatch(graphType, "All"))
			// Average Histograms graphs
			if(StringMatch(winName, "AvgDist_*"))  // Average Distance
				include = 1
			endif
			if(StringMatch(winName, "AvgInt_*"))  // Average Intensity
				include = 1
			endif
			if(StringMatch(winName, "AvgStep_*"))  // Average Stepsize
				include = 1
			endif
			if(StringMatch(winName, "AvgOntime_*"))  // Average On-time
				include = 1
			endif
			if(StringMatch(winName, "AvgOnrate_*"))  // Average On-rate
				include = 1
			endif
			if(StringMatch(winName, "Average*_C*E*"))  // Other average graphs
				include = 1
			endif
		endif
		
		if(StringMatch(graphType, "Compare") || StringMatch(graphType, "All"))
			// Compare Parameters graphs
			if(StringMatch(winName, "ColCmp_*"))  // Colocalization Compare
				include = 1
			endif
			if(StringMatch(winName, "Compare_*_C*E*"))  // Compare graphs with Col suffix
				include = 1
			endif
		endif
		
		if(include)
			graphList += winName + ";"
		endif
	endfor
	
	return graphList
End

// Colocalization Layout
Function CreateColocalizationLayout(pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	// Colocalization
	String graphList = GetColocalizationGraphList("All")
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs == 0)
		DoAlert 0, "No Colocalization graphs found.\nPlease run Colocalization analysis first."
		return 0
	endif
	
	Printf "Found %d Colocalization graphs\r", numGraphs
	
	// CreateAutoLayoutWithList
	return CreateAutoLayoutWithList("Colocalization", pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

// Batch Analysis 
Function CreateColBatchLayout(pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	String graphList = GetColocalizationGraphList("Batch")
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs == 0)
		DoAlert 0, "No Batch Analysis graphs found."
		return 0
	endif
	
	Printf "Found %d Batch Analysis graphs\r", numGraphs
	return CreateAutoLayoutWithList("ColBatch", pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

// Average Histograms 
Function CreateColAverageLayout(pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	String graphList = GetColocalizationGraphList("Average")
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs == 0)
		DoAlert 0, "No Average Histogram graphs found."
		return 0
	endif
	
	Printf "Found %d Average Histogram graphs\r", numGraphs
	return CreateAutoLayoutWithList("ColAverage", pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

// Compare Parameters 
Function CreateColCompareLayout(pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	String graphList = GetColocalizationGraphList("Compare")
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs == 0)
		DoAlert 0, "No Compare Parameters graphs found."
		return 0
	endif
	
	Printf "Found %d Compare Parameters graphs\r", numGraphs
	return CreateAutoLayoutWithList("ColCompare", pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

// =============================================================================
// Segmentation Layout Functions
// =============================================================================

// Segmentation
// graphType: "Batch", "Compare", "All"
Function/S GetSegmentationGraphList(graphType)
	String graphType
	
	String graphList = ""
	String allWindows = WinList("*", ";", "WIN:1")  // All graphs
	Variable numWindows = ItemsInList(allWindows)
	Variable i
	String winName
	
	for(i = 0; i < numWindows; i += 1)
		winName = StringFromList(i, allWindows)
		
		// 
		Variable include = 0
		
		if(StringMatch(graphType, "Batch") || StringMatch(graphType, "All"))
			// Batch Analysis graphs (_Seg* suffix)
			if(StringMatch(winName, "*_Seg*"))
				// Segment
				if(StringMatch(winName, "IntHist_*_Seg*"))  // Intensity histogram
					include = 1
				endif
				if(StringMatch(winName, "StepHist_*_Seg*"))  // Stepsize histogram
					include = 1
				endif
				if(StringMatch(winName, "MSD_*_Seg*"))  // MSD
					include = 1
				endif
				if(StringMatch(winName, "Duration_*_Seg*"))  // On-time
					include = 1
				endif
				if(StringMatch(winName, "Onrate_*_Seg*"))  // On-rate
					include = 1
				endif
			endif
		endif
		
		if(StringMatch(graphType, "Compare") || StringMatch(graphType, "All"))
			// Compare Parameters graphs
			if(StringMatch(winName, "SegCmp_*"))  // Segmentation Compare
				include = 1
			endif
		endif
		
		if(include)
			graphList += winName + ";"
		endif
	endfor
	
	return graphList
End

// Segmentation Layout
Function CreateSegmentationLayout(pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	// Segmentation
	String graphList = GetSegmentationGraphList("All")
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs == 0)
		DoAlert 0, "No Segmentation graphs found.\nPlease run Segmentation analysis first."
		return 0
	endif
	
	Printf "Found %d Segmentation graphs\r", numGraphs
	
	// CreateAutoLayoutWithList
	return CreateAutoLayoutWithList("Segmentation", pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

// Batch Analysis 
Function CreateSegBatchLayout(pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	String graphList = GetSegmentationGraphList("Batch")
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs == 0)
		DoAlert 0, "No Batch Analysis graphs found."
		return 0
	endif
	
	Printf "Found %d Batch Analysis graphs\r", numGraphs
	return CreateAutoLayoutWithList("SegBatch", pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

// Compare Parameters 
Function CreateSegCompareLayout(pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode)
	Variable pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode
	
	String graphList = GetSegmentationGraphList("Compare")
	
	Variable numGraphs = ItemsInList(graphList)
	if(numGraphs == 0)
		DoAlert 0, "No Compare Parameters graphs found."
		return 0
	endif
	
	Printf "Found %d Compare Parameters graphs\r", numGraphs
	return CreateAutoLayoutWithList("SegCompare", pageW_inch, pageH_inch, offset_mm, gap_mm, divW, divH, outputMode, graphList)
End

