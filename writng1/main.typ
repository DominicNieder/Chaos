#import "template.typ": thesis, placeholder-image

#show: thesis.with(title: "Chaos at the intersect of Quantum Mechanics and Classical Mechanics", author: "Dominic Nieder")

#align(center)[
  #text(size: 22pt, weight: "bold")[
    Chaos at the intersect of Quantum Mechanics and Classical Mechanics
  ]
  #v(1cm)
  #text(size: 14pt)[Dominic Nieder]
  #v(1cm)
  #text(size: 11pt)[September 2026]
]

#v(2cm)
#text[
  This thesis investigates...
]

#pagebreak()
#outline()
#pagebreak()

#include "chapters/01henonheiles.typ"

#bibliography("../literature/references.bib", style: "apa")