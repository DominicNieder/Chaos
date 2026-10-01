// template.typ
#let thesis(title: "", author: "", body) = {
  set document(title: title, author: author)
  set page(
    paper: "a4",
    margin: (left: 2.5cm, right: 2.5cm, top: 2.5cm, bottom: 2.5cm),
    numbering: "1",
  )
  set text(font: "Libertinus Serif", size: 11pt, lang: "en")  // or "New Computer Modern"
  set par(justify: true, leading: 0.75em)
  set heading(numbering: "1.1")
  show heading.where(level: 1): it => {
    pagebreak(weak: true)
    v(1em)
    it
    v(0.5em)
  }
  body
}

#let draft = sys.inputs.at("draft", default: "false") == "true"

#let placeholder-image(path, width: 100%) = {
  if draft {
    rect(width: width, height: 4cm, fill: luma(230), stroke: 0.5pt + gray)[
      #align(center + horizon)[
        #text(size: 8pt, fill: gray)[#path]
      ]
    ]
  } else {
    image(path, width: width)
  }
}