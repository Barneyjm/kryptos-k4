#!/bin/sh
# Project Gutenberg texts used for English/German n-gram scoring
cd "$(dirname "$0")"
for b in 1342 2701 1661 98 84; do curl -sL -o g$b.txt https://www.gutenberg.org/cache/epub/$b/pg$b.txt; done
for b in 2229 22367 7205 2407 6498; do curl -sL -o de$b.txt https://www.gutenberg.org/cache/epub/$b/pg$b.txt; done
ln -sf ../data/panel.txt ../data/split.txt ../data/tableau.txt .
