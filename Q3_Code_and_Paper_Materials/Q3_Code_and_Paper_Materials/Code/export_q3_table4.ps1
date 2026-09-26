param([string]$ResultsDirectory=(Join-Path $PSScriptRoot 'results'))
$ErrorActionPreference='Stop'
$rows=Import-Csv -LiteralPath (Join-Path $ResultsDirectory 'q3_final_energy_table.csv')
if($rows.Count -ne 2){throw 'Expected two verified strategy rows.'}
$caption='表4 不同辅助冷启动策略的优化结果及启动性能对比'
$headers=@('辅助冷启动策略','加热功率密度分配/W·cm⁻²','加热持续时间/s','各片辅助加热能耗/J','总辅助加热能耗/J','总启动时间/s','最大冰体积分数','启动结果')
$tableRows=@()
$body=''
$markdown="$caption`r`n`r`n|"+($headers -join '|')+"|`r`n|"+(('---','---','---:','---','---:','---:','---:','---') -join '|')+"|`r`n"
foreach($row in $rows){
    if($row.success -ne '1'){throw 'An unverified row cannot enter Table 4.'}
    $name=if($row.mode -eq 'preheat'){'纯预加热启动'}else{'恒定功率协同启动'}
    $q=1..5 | ForEach-Object {([double]$row."q$_").ToString('0.#########',[Globalization.CultureInfo]::InvariantCulture)}
    $energies=1..5 | ForEach-Object {([double]$row."E$($_)J").ToString('F3',[Globalization.CultureInfo]::InvariantCulture)}
    $values=@($name,('('+($q -join ', ')+')'),('{0:F4}' -f [double]$row.thUsedS),
        ('('+($energies -join ', ')+')'),('{0:F3}' -f [double]$row.ETotalJ),
        ('{0:F4}' -f [double]$row.startupTimeS),('{0:F6}' -f [double]$row.maximumIceVolumeFraction),'成功')
    $record=[ordered]@{}
    if($values.Count -ne 8){throw 'Table 4 must contain exactly eight columns.'}
    for($i=0;$i -lt 8;$i++){$record[$headers[$i]]=$values[$i]}
    $tableRows += [pscustomobject]$record
    $markdown+='|'+($values -join '|')+"|`r`n"
    $body+='<tr>'
    for($i=0;$i -lt 8;$i++){
        $value=[Net.WebUtility]::HtmlEncode($values[$i])
        if($i -eq 1){$value='('+($q[0..1] -join ', ')+',<br>'+($q[2..3] -join ', ')+',<br>'+$q[4]+')'}
        if($i -eq 3){$value='('+($energies[0..1] -join ', ')+',<br>'+($energies[2..3] -join ', ')+',<br>'+$energies[4]+')'}
        $body+='<td>'+$value+'</td>'
    }
    $body+='</tr>'
}
$note='注：五片共用一个加热持续时间。初温与环境温度均为−30 ℃，gammaIce=3.5，kFreeze=0.4，沿用原空间网格及0.01 K过零数值裕量。表内为当前搜索找到并复核通过的最低能耗候选，尚无全局最优性证明。各片能耗与总能耗分别四舍五入，完整精度见原始结果 CSV。'
$markdown+="`r`n$note`r`n"
$head='<tr>'+ (($headers | ForEach-Object {'<th>'+[Net.WebUtility]::HtmlEncode($_)+'</th>'}) -join '')+'</tr>'
$html=@"
<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>$caption</title>
<style>@page{size:A4 landscape;margin:18mm}body{font-family:SimSun,"Songti SC",serif;margin:32px;color:#000}table{border-collapse:collapse;table-layout:fixed;width:100%;max-width:1150px;margin:auto;font-size:12pt}caption{font-size:14pt;margin-bottom:12px}td,th{border:1px solid #000;text-align:center;vertical-align:middle;padding:10px 5px;line-height:1.65;overflow-wrap:anywhere}th{font-weight:bold}p{max-width:1150px;margin:16px auto;font-size:10pt;line-height:1.7}</style></head><body>
<table><caption>$caption</caption><colgroup><col style="width:11%"><col style="width:19%"><col style="width:9%"><col style="width:22%"><col style="width:10%"><col style="width:9%"><col style="width:11%"><col style="width:9%"></colgroup><thead>$head</thead><tbody>$body</tbody></table><p>$note</p></body></html>
"@
$stem=Join-Path $ResultsDirectory '表4_辅助冷启动优化结果'
$tableRows | Export-Csv -LiteralPath "$stem.csv" -NoTypeInformation -Encoding UTF8
[IO.File]::WriteAllText("$stem.md",$markdown,[Text.UTF8Encoding]::new($true))
[IO.File]::WriteAllText("$stem.html",$html,[Text.UTF8Encoding]::new($true))
Write-Output "Table 4 exported: $stem"
