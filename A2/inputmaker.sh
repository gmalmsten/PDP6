a=98


rm -f input$a.txt
echo -n "$a " >> input$a.txt
for i in $(seq 1 $((2*a*a)))
do
    echo -n "$i " >> input$a.txt
done
