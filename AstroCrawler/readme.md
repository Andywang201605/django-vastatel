### Astronomical Notice Parser
This directory contains different (only Atel so far) crawlers for astronomical notice

- `coordMatch.py`: Different coordinate matcher/searcher and coordinate filter (i.e., give unique)
    - `coordMatch`: 
        Use regular expression to match different format coordinates in a given text. \
        `>>> coordmatcher = CoordMatch(rawString)`\
        `>>> coordStrings, coordRaws = coordMatch.matchAll()`\
        `coordStrings` is a list contains all parsed coordinates, you can convert them to astropy objects simply with `SkyCoord(coordStrings)`; `coordRaws` is a list contains all matched text, mainly used as a reference.
    - `TNSQuery`:
        For a given Transient Name Server (TNS) object name, it will use `requests` to fetch coordinate information. Note: we use normal requests to get information from TNS, but it is recommended to use API to query information. \
        `>>> tnsquery = TNSQuery('2021yte')` \
        `>>> coordStr = tnsquery.checkcoord()` \
        `coordStr` is a string with a format `HH:MM:SS.SSS DD:MM:SS.SS`, you can convert it to an astropy object with `SkyCoord(coordStr, unit=(u.hourangle, u.degree))`
    - `TitleSearch`:
        It will try searching any object names in the title (or a short sentence). It will split the whole sentence into different fractions and use Simbad to search for coordinates. With `withPrefix=True`(by default `False`), it will only looking for words after certain words: of, from, in, blazer, star, and nova. Note: `withPrefix=True` will miss object at the beginning of the sentence. \
        `>>> titlesearch = TitleSearch(title, phraseMaxWord=3)` \
        `>>> titleCoords, titleNames = titlesearch.simbadsearch(withPrefix=False, removeThreshold=2)` \
        `titleCoords` is a list contains all parsed coordinates, you can convert them to astropy objects with `SkyCoord(coordStr, unit=(u.hourangle, u.degree))`; `titleNames` is a list contains all matched object name. \
        Note: Some objects may have weired names (for example: A1), there are chances that it will find something irrelevant; `removeThreshold` is used for remove some names with lots of astronomical objects matched (e.g., MAGIC).
    - `CoordFilter`:
        Given a list of coordinates (either a `SkyCoord` object, or a list of `SkyCoord` objects), it will remove all coordinates that are too close and leave a single one. \
        `>>> coordfilter = CoordFilter(coordLst)` \
        `>>> coordfilter.filterSource(radius=5.)` \
        You can access the final result by `coordfilter.uniqueSources` - it is a list with coordinates in decimal places, you can convert them to astropy objects with `SkyCoord(coordfilter.uniqueSources, unit=u.degree)`

<hr />

#### ATel
Scripts for monitor, download, parse and send message for Astronomer Telegrams (https://astronomerstelegram.org/)

