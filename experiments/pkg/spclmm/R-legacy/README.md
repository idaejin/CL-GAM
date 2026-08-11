# Quarantined from R/ in experiment fork 0.1.1
# These files call retired rgeos (gTouches, area.poly, intersect).
# Paper Madrid analyses use precomputed C matrices (spC_*.Rdata) + pois_SAP;
# they do not need atp_setup*.
# TODO: rewrite with sf::st_intersects / st_area if ATP construction is revived.
