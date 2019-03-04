library(DiagrammeR)
library(DiagrammeRsvg)
library(V8)
library(convertGraph)
install.phantom("/home/hoanguc3m/bin/phantomjs")

graph <- grViz("
               digraph subscript {

               # Graph statements
               graph [layout = dot, overlap = false]

               # Node statements
               node [shape = circle,
               style = filled,
               color = black,
               penwidth = 2.0,
               fillcolor = orange]

               node [fillcolor = OrangeRed]
               a [label =<V<SUB>0</SUB>>];

               node [shape = circle,
               penwidth = 2.0,
               fillcolor = YellowGreen]
               b1 [label =<U<SUB>1</SUB>>];
               b2 [label =<U<SUB>2</SUB>>];
               b3 [label =<U<SUB>3</SUB>>];
               b4 [label =<U<SUB>4</SUB>>];
               b5 [label =<U<SUB>5</SUB>>];
               # Edge statements
               edge [color = black, arrowhead = none]
               a -> {b1, b2, b3, b4, b5}
               }
               ")

cat(export_svg(graph), file='onefactor.svg')
convertGraph("onefactor.svg", "onefactor.pdf")

graph <- grViz("
               digraph subscript {

               # Graph statements
               graph [layout = dot, overlap = false]

               # Node statements
               node [shape = circle,
               style = filled,
               color = black,
               penwidth = 2.0,
               fillcolor = orange]

               node [fillcolor = OrangeRed]
               a1 [label =<V<SUB>0</SUB>>];
               a2 [label =<V<SUB>1</SUB>>];

               subgraph { rank = same; a1}

               node [shape = circle,
               penwidth = 2.0,
               fillcolor = YellowGreen]
               b1 [label =<U<SUB>1</SUB>>];
               b2 [label =<U<SUB>2</SUB>>];
               b3 [label =<U<SUB>3</SUB>>];
               b4 [label =<U<SUB>4</SUB>>];
               b5 [label =<U<SUB>5</SUB>>];
               # Edge statements
               edge [color = black, arrowhead = none]
               a1 -> {b1, b2, b3, b4, b5}
               a2 -> {b1, b2, b3, b4, b5}
               edge [color = white, arrowhead = none]
               a1 -> a2
               }
               ")
cat(export_svg(graph), file='twofactors.svg')
convertGraph("twofactors.svg", "twofactors.pdf")


graph <- grViz("
               digraph subscript {

               # Graph statements
               graph [layout = dot, overlap = false]

               # Node statements
               node [shape = circle,
               style = filled,
               color = black,
               penwidth = 2.0,
               fillcolor = orange]

               node [fillcolor = OrangeRed]
               a0 [label =<V<SUB>0</SUB>>];
               a1 [label =<V<SUB>1</SUB>>];
               a2 [label =<V<SUB>2</SUB>>];
               a3 [label =<V<SUB>3</SUB>>];

               node [shape = circle,
               penwidth = 2.0,
               fillcolor = YellowGreen]
               b1 [label =<U<SUB>1</SUB>>];
               b2 [label =<U<SUB>2</SUB>>];
               b3 [label =<U<SUB>3</SUB>>];
               b4 [label =<U<SUB>4</SUB>>];
               b5 [label =<U<SUB>5</SUB>>];
               b6 [label =<U<SUB>6</SUB>>];
               b7 [label =<U<SUB>7</SUB>>];
               b8 [label =<U<SUB>8</SUB>>];
               b9 [label =<U<SUB>9</SUB>>];
               b10 [label =<U<SUB>10</SUB>>];
               b11 [label =<U<SUB>11</SUB>>];
               b12 [label =<U<SUB>12</SUB>>];
               # Edge statements
               edge [color = white, arrowhead = none]
               a0 -> {a1, a2, a3}

               edge [color = black, arrowhead = none]
               a0 -> {b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12}
               a1 -> {b1, b2, b3, b4}
               a2 -> {b5, b6, b7, b8}
               a3 -> {b9, b10, b11, b12}
               }
               ")
graph
cat(export_svg(graph), file='bifactors.svg')
convertGraph("bifactors.svg", "bifactors.pdf")

graph <- grViz("
               digraph subscript {

               # Graph statements
               graph [layout = dot, overlap = false]

               # Node statements
               node [shape = circle,
               style = filled,
               color = black,
               penwidth = 2.0,
               fillcolor = orange]

               node [fillcolor = OrangeRed]
               a0 [label =<V<SUB>0</SUB>>];
               a1 [label =<V<SUB>1</SUB>>];
               a2 [label =<V<SUB>2</SUB>>];
               a3 [label =<V<SUB>3</SUB>>];

               node [shape = circle,
               penwidth = 2.0,
               fillcolor = YellowGreen]
               b1 [label =<U<SUB>1</SUB>>];
               b2 [label =<U<SUB>2</SUB>>];
               b3 [label =<U<SUB>3</SUB>>];
               b4 [label =<U<SUB>4</SUB>>];
               b5 [label =<U<SUB>5</SUB>>];
               b6 [label =<U<SUB>6</SUB>>];
               b7 [label =<U<SUB>7</SUB>>];
               b8 [label =<U<SUB>8</SUB>>];
               b9 [label =<U<SUB>9</SUB>>];
               b10 [label =<U<SUB>10</SUB>>];
               b11 [label =<U<SUB>11</SUB>>];
               b12 [label =<U<SUB>12</SUB>>];
               # Edge statements
               edge [color = black, arrowhead = none]
               a0 -> {a1, a2, a3}
               a1 -> {b1, b2, b3, b4}
               a2 -> {b5, b6, b7, b8}
               a3 -> {b9, b10, b11, b12}
               }
               ")
cat(export_svg(graph), file='nestedfactors.svg')
convertGraph("nestedfactors.svg", "nestedfactors.pdf")



graph <- grViz("
               digraph subscript {

               # Graph statements
               graph [layout = dot, overlap = false, rankdir = TD]

               # Node statements

               node [shape = circle,
               style = filled,
               color = black,
               penwidth = 2.0,
               fillcolor = YellowGreen]
               b1 [label =<1>];
               b2 [label =<2>];
               b3 [label =<3>];
               b4 [label =<4>];

               b5 [label =<12>];
               b6 [label =<23>];
               b7 [label =<34>];

               b8 [label =<13|2>];
               b9 [label =<24|3>];


               node [shape = plaintext,
               style = none,
               color = black,
               penwidth = 2.0]
               T1 [label =<T<SUB>1</SUB>>];
               T2 [label =<T<SUB>2</SUB>>];
               T3 [label =<T<SUB>3</SUB>>];
               # {rank = same; T1; T2; T3}

               # Edge statements
               edge [color = white, arrowhead = none]
               T1 -> b1
               #    T1 -> T2
               T2 -> b5
               #    T2 -> T3
               T3 -> b8

               edge [color = black, arrowhead = none]
               b1 -> b2 [label =<12>]
               b2 -> b3 [label =<23>]
               b3 -> b4 [label =<34>]
               b5 -> b6 [label =<13|2>]
               b6 -> b7 [label =<24|3>]
               b8 -> b9 [label =<14|23>]
               }
               ")
graph
cat(export_svg(graph), file='Dvine.svg')
convertGraph("Dvine.svg", "Dvine.pdf")



graph <- grViz("
               digraph subscript {

               # Graph statements
               graph [layout = dot, overlap = false, rankdir = TD]

               # Node statements

               node [shape = circle,
               style = filled,
               color = black,
               penwidth = 2.0,
               fillcolor = YellowGreen]
               b1 [label =<1>];
               b2 [label =<2>];
               b3 [label =<3>];
               b4 [label =<4>];

               b5 [label =<12>];
               b6 [label =<13>];
               b7 [label =<14>];

               b8 [label =<23|1>];
               b9 [label =<24|1>];


               node [shape = plaintext,
               style = none,
               color = black,
               penwidth = 2.0]
               T1 [label =<T<SUB>1</SUB>>];
               T2 [label =<T<SUB>2</SUB>>];
               T3 [label =<T<SUB>3</SUB>>];


               # Edge statements
               edge [color = white, arrowhead = none]
               T1 -> b1
               #    T1 -> T2
               T2 -> b5
               #    T2 -> T3
               T3 -> b8

               edge [color = black, arrowhead = none]
               b1 -> b2 [label =<12>]
               b1 -> b3 [label =<13>]
               b1 -> b4 [label =<14>]
               b5 -> b6 [label =<23|1>]
               b5 -> b7 [label =<24|1>]
               b8 -> b9 [label =<34|12>]
               }
               ")
graph
cat(export_svg(graph), file='Cvine.svg')
convertGraph("Cvine.svg", "Cvine.pdf")

set.seed(32323)
graph <- grViz("
    digraph subscript {

               # Graph statements
               graph [layout = dot, overlap = false, splines = true]

               # Node statements

               node [shape = circle,
               style = filled,
               color = black,
               penwidth = 2.0,
               fillcolor = YellowGreen]
               b1 [label =<1>];
               b2 [label =<2>];
               b3 [label =<3>];
               b4 [label =<4>];

               subgraph { rank = same; b2; b3;}

               edge [color = black, arrowhead = none]
               b1 -> b2 [label =<12>]
               b1 -> b3 [label =<13>]
               b1 -> b4 [label =<14>]
               b2 -> b3 [label =<23|1>]
               b2 -> b4 [label =<24|1>]
               b3 -> b4 [label =<34|12>]
               }
               ")
graph
cat(export_svg(graph), file='Cvine_full.svg')
convertGraph("Cvine_full.svg", "Cvine_full.pdf")
