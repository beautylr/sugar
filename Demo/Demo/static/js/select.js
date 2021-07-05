var selects;
    window.onload = function () {
        var form = document.getElementById("i_form");
        selects = form.getElementsByTagName("select");
        for (var i = 0; i < 5; i++) {
            selects[i].addEventListener("change",function(){change(this);});
        }

    }

    function change(e){
        var t = true;
        var id = parseInt(e.id);
        for (var i = 0; i < 5; i++){
            if(t){
                if(id == i){
                    t = false;
                }
            }else{
                selects[i].innerHTML = "<option>0.0</option>";
                var j = "0.1";
                var gross = updateGross();
                while ( j <= gross ){
                    var op = document.createElement("option");
                    op.innerText = j;
                    selects[i].appendChild(op);
                    j = addtion(j,"0.1");
                }

            }
        }
    }

    function  updateGross(){
        var sum = 0;
        for (var x = 0; x < 5; x++){
            var index = selects[x].selectedIndex;
            var text = selects[x].options[index].text;
            var value = parseFloat(text);
            sum += value;
        }
        var gross = 1.0 - sum;
        return gross.toFixed(1);
    }

    function addtion(a, b) {
        var sum = parseFloat(a)+parseFloat(b);
        return sum.toFixed(1);
    }