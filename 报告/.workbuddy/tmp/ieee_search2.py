import urllib.request, urllib.parse, json, time

KEY = "4sbyaxjws384aw362s3t5txx"
BASE = "https://ieeexploreapi.ieee.org/api/v1/search/articles"

queries = [
    'PEMFC cold start optimization',
    '"fuel cell stack" "cold start"',
    'fuel cell subfreezing startup',
    'PEM fuel cell cold start preheating energy',
    'fuel cell cold start experiment voltage temperature',
    'PEMFC freeze start model validation',
]

for q in queries:
    params = {"apikey": KEY, "format": "json", "querytext": q, "max_records": 8}
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
        print(f"  摘要: {ab[:450]}")
    time.sleep(1.2)
