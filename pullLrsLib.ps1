$url = "http://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/"

$html = curl -s $url


$htmlPattern = '>lrslib-([^\.]+)\.tar\.gz[^\d]*(\d{4}\-\d{2}\-\d{2}\s+\d{2}:\d{2})'

$latestDate = [System.DateTime]"2000-01-01"
$latestVersion = ""

foreach ($line in ($html.Split("\n"))) {
    $match = [Regex]::Match($line,$htmlPattern)
    
    if ($match.Success) {
        $currentDate = [System.DateTime]$match.Groups[2].Value
        if ($currentDate -gt $latestDate) {
            $latestVersion = $match.Groups[1].Value
            $latestDate = $currentDate
        }
    }
    
}

$baseVersion = [Regex]::match($latestVersion,"\d+").Value

curl --output "lrslib.tar.gz" ($url + "lrslib-" + $latestVersion + ".tar.gz")

tar xvzf ./lrslib.tar.gz