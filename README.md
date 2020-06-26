# p4ppk
Phantom 4 RTK PPK processing tools

Several tools to help with adding PPK trajectory solutions in to Phantom 4 RTK surveys. Meant to work with RTKLIB .pos trajectory solutions and 

avgPos.py - Calculate average position from a .pos file.

ppkMerge.py - Ingest the DJI timestamp file and PPK trajectory (.pos) and calculate new positions for each photo - output a CSV with the new position.

writeExifXmp.py - Write corrected positions to the XMP and EXIF data in the survey images.
