param([string]$ResultsDirectory=(Join-Path $PSScriptRoot 'results'))
$ErrorActionPreference='Stop'
Add-Type -AssemblyName System.Drawing
$stem=Join-Path $ResultsDirectory '表4_辅助冷启动优化结果'
$word=$null; $document=$null
try {
    $word=New-Object -ComObject Word.Application
    $word.Visible=$false
    $word.DisplayAlerts=0
    $document=$word.Documents.Add()
    $document.Content.InsertAfter("表4 不同辅助冷启动策略的优化结果及启动性能对比`r")
    $document.PageSetup.PaperSize=7
    $document.PageSetup.Orientation=1
    $document.PageSetup.TopMargin=36
    $document.PageSetup.BottomMargin=36
    $document.PageSetup.LeftMargin=36
    $document.PageSetup.RightMargin=36
    $document.Content.Font.NameFarEast='宋体'
    $document.Content.Font.Name='Times New Roman'
    $document.Content.Font.Size=11
    $data=Import-Csv -LiteralPath "$stem.csv" -Encoding UTF8
    $headers=@($data[0].PSObject.Properties.Name)
    if($headers.Count -ne 8 -or $data.Count -ne 2){throw 'Expected an eight-column, two-strategy table.'}
    $range=$document.Range($document.Content.End-1,$document.Content.End-1)
    $table=$document.Tables.Add($range,3,8)
    for($i=1;$i -le 8;$i++){
        $table.Cell(1,$i).Range.Text=$headers[$i-1]
        for($r=2;$r -le 3;$r++){
            $value=[string]$data[$r-2].($headers[$i-1])
            if($i -eq 2 -or $i -eq 4){
                $parts=$value.Trim([char[]]'()').Split(',') | ForEach-Object {$_.Trim()}
                $value='('+($parts[0..1] -join ', ')+','+[char]11+($parts[2..3] -join ', ')+','+[char]11+$parts[4]+')'
            }
            $table.Cell($r,$i).Range.Text=$value
        }
    }
    $table.AllowAutoFit=$false
    $width=$document.PageSetup.PageWidth-$document.PageSetup.LeftMargin-$document.PageSetup.RightMargin
    $fractions=@(.11,.19,.09,.22,.10,.09,.11,.09)
    for($r=1;$r -le $table.Rows.Count;$r++){
        for($i=1;$i -le 8;$i++){$table.Cell($r,$i).SetWidth($width*$fractions[$i-1],0)}
    }
    $table.Rows.AllowBreakAcrossPages=$false
    $table.Range.ParagraphFormat.Alignment=1
    $table.Range.Cells.VerticalAlignment=1
    $table.TopPadding=8
    $table.BottomPadding=8
    $table.Borders.Enable=1
    $table.Borders.OutsideColor=0
    $table.Borders.InsideColor=0
    foreach($border in @(-1,-2,-3,-4,-5,-6)){
        $table.Borders.Item($border).LineStyle=1
        $table.Borders.Item($border).LineWidth=4
        $table.Borders.Item($border).Color=0
    }
    $table.Rows.Item(1).HeadingFormat=-1
    $document.Content.InsertAfter("`r注：五片共用一个加热持续时间；初温和环境温度均为−30 ℃，gammaIce=3.5，kFreeze=0.4。沿用原空间网格和0.01 K过零数值裕量。表内为当前已复核的最低能耗可行候选，尚无全局最优性证明。各项能耗分别四舍五入，完整精度见原始结果CSV。")
    $document.Content.Font.Bold=0
    for($i=1;$i -le 8;$i++){$table.Cell(1,$i).Range.Font.Bold=-1}
    $document.Paragraphs.Item(1).Range.Font.Size=14
    $document.Paragraphs.Item(1).Range.ParagraphFormat.Alignment=1
    $document.Repaginate()
    $docx="$stem.docx"; $format=16
    $document.SaveAs2([ref]$docx,[ref]$format)
    $document.ExportAsFixedFormat("$stem.pdf",17)
    $document.ActiveWindow.View.Type=3
    $pages=$document.ActiveWindow.Panes.Item(1).Pages
    if($pages.Count -ne 1){throw "Expected one table page, got $($pages.Count)."}
    # Windows PDF rendering preserves the printed borders and font weights.
    $renderer=Join-Path $PSScriptRoot 'render_q3_table4_pdf.ps1'
    & powershell -NoProfile -ExecutionPolicy Bypass -File $renderer -ResultsDirectory $ResultsDirectory
    if($LASTEXITCODE -ne 0){throw 'PDF visual verification render failed.'}
    Copy-Item -LiteralPath (Join-Path $ResultsDirectory 'table4_pdf_render.png') -Destination (Join-Path $ResultsDirectory 'table4_word_preview.png') -Force
    Write-Output "Editable Table 4 and one-page PDF exported: $docx"
} finally {
    if($null -ne $document){$document.Close(0)}
    if($null -ne $word){if($word.Documents.Count -eq 0){$word.Quit()};[void][Runtime.InteropServices.Marshal]::ReleaseComObject($word)}
}
