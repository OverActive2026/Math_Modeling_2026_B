import urllib.request, urllib.parse, json, time

KEY = "4sbyaxjws384aw362s3t5txx"
BASE = "https://ieeexploreapi.ieee.org/api/v1/search/articles"

queries = [
    'fuel cell purge water removal',
    'PEM fuel cell purge shutdown',
    '"fuel cell" "cold start" purge',
    'fuel cell water evacuation cold start',
    'fuel cell nitrogen purge drying',
    'PEMFC purge strategy residual water',
    'fuel cell cathode purge ice',
    'fuel cell shutdown purge experiment',
]

for q in queries:
    params = {"apikey": KEY, "format": "json", "querytext": q, "max_records": 10}
    url = BASE + "?" + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.load(r)
    except Exception as e:
        print(f"### QUERY: {q}  -> ERROR {e}")
        time.sleep(1.2)
        continue
    arts = data.get("articles", [])
    print(f"\n### QUERY: {q}  (total={data.get('total_records')}, shown={len(arts)})")
    for a in arts:
        authors = ", ".join(x.get("full_name", "") for x in a.get("authors", {}).get("authors", [])[:3])
        print(f"- [{a.get('publication_year')}] {a.get('title')}")
        print(f"  期刊: {a.get('publication_title')} | 作者: {authors} | DOI: {a.get('doi')}")
        ab = (a.get("abstract") or "").replace("\n", " ")
        print(f"  摘要: {ab[:500]}")
    time.sleep(1.2)
