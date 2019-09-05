pandoc --from markdown --to latex -o user_guide.pdf README.md --template=template.tex
echo '# CIAlign' > ../README.md
awk 'NR > 7' README.md >> ../README.md
