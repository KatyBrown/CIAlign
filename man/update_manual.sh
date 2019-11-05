pandoc --from markdown --to latex -o user_guide.pdf user_guide.md --template=template.tex
echo '# CIAlign' > ../README.md
awk 'NR > 7' user_guide.md >> ../README.md
