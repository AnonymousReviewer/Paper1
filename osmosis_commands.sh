
# besancon area extraction -- using pbf file from geofabrik downloaded on Mar 13 sep 2016 16:46:31 CEST 
osmosis --read-pbf data/franche-comte-latest.osm.pbf --log-progress --bounding-box left=5.9191 right=6.0639 bottom=47.2038 top=47.281 --tf accept-ways building=* --used-node --write-xml buildings_besac.osm


