param([string]$ResultsDirectory=(Join-Path $PSScriptRoot 'results'))
$ErrorActionPreference='Stop'
Add-Type -AssemblyName System.Runtime.WindowsRuntime
[Windows.Storage.StorageFile,Windows.Storage,ContentType=WindowsRuntime] | Out-Null
[Windows.Data.Pdf.PdfDocument,Windows.Data.Pdf,ContentType=WindowsRuntime] | Out-Null
[Windows.Data.Pdf.PdfPageRenderOptions,Windows.Data.Pdf,ContentType=WindowsRuntime] | Out-Null
[Windows.Storage.Streams.InMemoryRandomAccessStream,Windows.Storage.Streams,ContentType=WindowsRuntime] | Out-Null
function Await-Result($Operation,[Type]$ResultType){
    $method=[System.WindowsRuntimeSystemExtensions].GetMethods() | Where-Object {
        $_.Name -eq 'AsTask' -and $_.IsGenericMethod -and $_.GetGenericArguments().Count -eq 1 -and $_.GetParameters().Count -eq 1
    } | Select-Object -First 1
    $task=$method.MakeGenericMethod($ResultType).Invoke($null,@($Operation))
    $task.Wait();return $task.Result
}
$path=(Resolve-Path -LiteralPath (Join-Path $ResultsDirectory '表4_辅助冷启动优化结果.pdf')).Path
$file=Await-Result ([Windows.Storage.StorageFile]::GetFileFromPathAsync($path)) ([Windows.Storage.StorageFile])
$pdf=Await-Result ([Windows.Data.Pdf.PdfDocument]::LoadFromFileAsync($file)) ([Windows.Data.Pdf.PdfDocument])
if($pdf.PageCount -ne 1){throw "Expected one page, got $($pdf.PageCount)."}
$page=$pdf.GetPage(0)
$stream=[Windows.Storage.Streams.InMemoryRandomAccessStream]::new()
$options=[Windows.Data.Pdf.PdfPageRenderOptions]::new();$options.DestinationWidth=1600
$operation=$page.RenderToStreamAsync($stream,$options)
$actionMethod=[System.WindowsRuntimeSystemExtensions].GetMethods() | Where-Object {
    $_.Name -eq 'AsTask' -and -not $_.IsGenericMethod -and $_.GetParameters().Count -eq 1 -and $_.GetParameters()[0].ParameterType.Name -eq 'IAsyncAction'
} | Select-Object -First 1
$task=$actionMethod.Invoke($null,@($operation));$task.Wait()
$read=[System.IO.WindowsRuntimeStreamExtensions]::AsStreamForRead($stream.GetInputStreamAt(0))
$output=[IO.File]::Create((Join-Path $ResultsDirectory 'table4_pdf_render.png'))
$read.CopyTo($output);$output.Dispose();$read.Dispose();$stream.Dispose();$page.Dispose()
Write-Output 'Word-exported PDF rendered and page count checked.'
