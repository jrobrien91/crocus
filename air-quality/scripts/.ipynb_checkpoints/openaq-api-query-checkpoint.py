import requests

chicago_bbox = [-88.885961, 41.274839, -87.117162, 42.477288]


# Check the number of locations within the "Chicago-Naperville-Joliet" city
##url = "https://api.openaq.org/v2/cities?limit=1000&page=1&offset=0&sort=asc&country_id=13&country=US&city=Chicago-Naperville-Joliet&order_by=city"

##headers = {"accept": "application/json"}

##response = requests.get(url, headers=headers)

##print(response.text)

# Get Locations and Observations

url = "https://api.openaq.org/v2/locations?limit=1000&page=1&offset=0&sort=desc&parameter_id=2&radius=1000&country=US&city=Chicago-Naperville-Joliet&order_by=lastUpdated&dump_raw=false"

headers = {"accept": "application/json"}

response = requests.get(url, headers=headers)

print(response.text)
