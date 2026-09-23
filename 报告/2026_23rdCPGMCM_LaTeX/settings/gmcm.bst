ENTRY
  { address
    author
    booktitle
    editor
    howpublished
    institution
    journal
    language
    note
    number
    organization
    pages
    publisher
    school
    title
    urldate
    url
    volume
    year
  }
  {}
  {}

INTEGERS { nameptr numnames totalentries }
STRINGS { names result item sep colon period etal }

FUNCTION {is.english}
{ language empty$
    { #0 }
    { language "en" =
        { #1 }
        { language "english" = }
      if$
    }
  if$
}

FUNCTION {set.punctuation}
{ is.english
    { ", " 'sep :=  ": " 'colon :=  "." 'period :=  "et al." 'etal := }
    { "，" 'sep :=  "：" 'colon :=  "。" 'period :=  "等" 'etal := }
  if$
}

FUNCTION {field.or.empty}
{ duplicate$ empty$ { pop$ "" } 'skip$ if$ }

FUNCTION {append}
{ 'item :=
  item empty$
    'skip$
    { result empty$
        { item 'result := }
        { result sep * item * 'result := }
      if$
    }
  if$
}

FUNCTION {format.names}
{ 'names :=
  #1 'nameptr :=
  names num.names$ 'numnames :=
  ""
  { nameptr numnames #1 + < }
  { names nameptr "{ff~}{vv~}{ll}{, jj}" format.name$ 'item :=
    item "others" =
      { etal 'item := }
      'skip$
    if$
    nameptr #1 >
      { sep * }
      'skip$
    if$
    item *
    nameptr #1 + 'nameptr :=
  }
  while$
}

FUNCTION {format.author}
{ author empty$
    { editor empty$
        { organization field.or.empty }
        { editor format.names }
      if$
    }
    { author format.names }
  if$
}

FUNCTION {start.entry}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  newline$
  "" 'result :=
  set.punctuation
  format.author append
  title field.or.empty append
  title empty$ { "Missing title in " cite$ * warning$ } 'skip$ if$
}

FUNCTION {finish.entry}
{ note field.or.empty append
  result period * write$
  newline$
}

FUNCTION {format.imprint}
{ address field.or.empty
  publisher empty$
    'skip$
    { duplicate$ empty$
        { pop$ publisher }
        { colon * publisher * }
      if$
    }
  if$
}

FUNCTION {format.volume.pages}
{ volume field.or.empty
  number empty$
    'skip$
    { "(" * number * ")" * }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ pages }
        { colon * pages * }
      if$
    }
  if$
}

FUNCTION {format.url}
{ url empty$
    { howpublished field.or.empty }
    { "\url{" url * "}" * }
  if$
}

FUNCTION {format.access}
{ urldate empty$
    { url empty$
        'skip$
        { "Missing urldate in " cite$ * warning$ }
      if$
      ""
    }
    { language empty$
        { "访问时间：" urldate * }
        { language "en" = language "english" = or
            { "Accessed: " urldate * }
            { "访问时间：" urldate * }
          if$
        }
      if$
    }
  if$
}

FUNCTION {book}
{ start.entry
  format.imprint append
  pages field.or.empty append
  year field.or.empty append
  finish.entry
}

FUNCTION {article}
{ start.entry
  journal field.or.empty append
  format.volume.pages append
  year field.or.empty append
  finish.entry
}

FUNCTION {inproceedings}
{ start.entry
  booktitle field.or.empty append
  format.imprint append
  pages field.or.empty append
  year field.or.empty append
  finish.entry
}

FUNCTION {incollection} { inproceedings }
FUNCTION {conference} { inproceedings }

FUNCTION {techreport}
{ start.entry
  institution field.or.empty append
  number field.or.empty append
  year field.or.empty append
  finish.entry
}

FUNCTION {phdthesis}
{ start.entry
  school field.or.empty append
  year field.or.empty append
  finish.entry
}

FUNCTION {mastersthesis} { phdthesis }

FUNCTION {misc}
{ start.entry
  format.url append
  url empty$
    { year field.or.empty append }
    { format.access append }
  if$
  finish.entry
}

FUNCTION {online} { misc }
FUNCTION {electronic} { misc }
FUNCTION {manual} { misc }
FUNCTION {unpublished} { misc }
FUNCTION {default.type} { misc }

FUNCTION {init.count}
{ #0 'totalentries := }

FUNCTION {count.entry}
{ totalentries #1 + 'totalentries := }

FUNCTION {begin.bib}
{ "\begin{thebibliography}{" write$
  totalentries int.to.str$ write$
  "}" write$
  newline$
}

FUNCTION {end.bib}
{ newline$ "\end{thebibliography}" write$ newline$ }

READ
EXECUTE {init.count}
ITERATE {count.entry}
EXECUTE {begin.bib}
ITERATE {call.type$}
EXECUTE {end.bib}
