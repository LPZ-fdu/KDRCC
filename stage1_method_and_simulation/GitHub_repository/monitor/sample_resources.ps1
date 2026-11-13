param([Parameter(Mandatory=$true)][string]$Folder,
      [Parameter(Mandatory=$true)][int]$OwnerPid,
      [double]$SampleSeconds=1,[double]$DiskSeconds=5)
# Process-local execution policy only; this script reads process statistics and
# task-owned context/cache paths. It does not change OS settings or user files.
$ErrorActionPreference='Stop'
$culture=[Globalization.CultureInfo]::InvariantCulture
$utf8=New-Object Text.UTF8Encoding($false)
$processWriter=$null
$totalWriter=$null
function CsvCell($value) {
  if($null -eq $value){return ''}
  if($value -is [IFormattable]){$value=$value.ToString($null,$culture)}
  '"'+([string]$value).Replace('"','""')+'"'
}
function WriteCsvLine($writer,$values){$writer.WriteLine((($values | ForEach-Object {CsvCell $_}) -join ','))}
function DirectoryBytes([string]$path){
  if([string]::IsNullOrEmpty($path) -or !(Test-Path -LiteralPath $path -PathType Container)){return 0.0}
  $sum=0.0
  Get-ChildItem -LiteralPath $path -File -Recurse -ErrorAction SilentlyContinue | ForEach-Object {$sum+=$_.Length}
  return $sum
}
try {
  $parent=Get-Process -Id $OwnerPid
  $ownerStart=$parent.StartTime
  $processWriter=New-Object IO.StreamWriter((Join-Path $Folder 'process_samples.csv'),$false,$utf8)
  $totalWriter=New-Object IO.StreamWriter((Join-Path $Folder 'total_samples.csv'),$false,$utf8)
  $processWriter.AutoFlush=$true;$totalWriter.AutoFlush=$true
  $processWriter.WriteLine('time_epoch,pid,role,kind,key,phase,N,L,M,kappa,lambda,working_set_bytes,private_bytes,cache_bytes')
  $totalWriter.WriteLine('time_epoch,process_count,summed_working_set_bytes,summed_private_bytes,total_cache_bytes')
  [IO.File]::WriteAllText((Join-Path $Folder 'ready.txt'),[string]$PID,$utf8)
  $lastDisk=[DateTime]::MinValue
  while(!(Test-Path -LiteralPath (Join-Path $Folder 'stop.txt'))){
    $owner=Get-Process -Id $OwnerPid -ErrorAction SilentlyContinue
    if($null -eq $owner -or $owner.StartTime -ne $ownerStart){break}
    $now=[DateTime]::UtcNow;$epoch=($now-[DateTime]'1970-01-01').TotalSeconds
    $sampleDisk=($now-$lastDisk).TotalSeconds -ge $DiskSeconds
    $rss=0.0;$private=0.0;$count=0;$totalDisk=$null
    $files=Get-ChildItem -LiteralPath (Join-Path $Folder 'contexts') -Filter '*.csv' -File
    foreach($file in $files){
      try {
        $c=Import-Csv -LiteralPath $file.FullName -Encoding UTF8 | Select-Object -First 1
        if($null -eq $c){continue}
        $p=Get-Process -Id ([int]$c.pid) -ErrorAction SilentlyContinue
        if($null -eq $p){continue}
        $p.Refresh();$w=[double]$p.WorkingSet64;$v=[double]$p.PrivateMemorySize64
        $cache=$null
        if($sampleDisk){$cache=DirectoryBytes $c.work;if($c.role -eq 'main'){$totalDisk=$cache}}
        WriteCsvLine $processWriter @($epoch,$c.pid,$c.role,$c.kind,$c.key,$c.phase,$c.N,$c.L,$c.M,$c.kappa,$c.lambda,$w,$v,$cache)
        $rss+=$w;$private+=$v;$count++
      } catch { }
    }
    WriteCsvLine $totalWriter @($epoch,$count,$rss,$private,$totalDisk)
    if($sampleDisk){$lastDisk=$now}
    Start-Sleep -Milliseconds ([int][Math]::Max(100,$SampleSeconds*1000))
  }
} catch {
  [IO.File]::WriteAllText((Join-Path $Folder 'error.txt'),($_ | Out-String),$utf8)
} finally {
  if($null -ne $processWriter){$processWriter.Close()}
  if($null -ne $totalWriter){$totalWriter.Close()}
  [IO.File]::WriteAllText((Join-Path $Folder 'finished.txt'),[DateTime]::UtcNow.ToString('o'),$utf8)
}
