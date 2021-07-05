
var dotDiv = document.getElementById("dot-chart");



var traceA = {

  type: "scatter",

  mode: "markers",

  x: [2011, 2012, 2013, 2014, 2015, 2016],

  y: [789, 795, 760, 775, 780, 783],

  name: 'Highest Marks',

  marker: {

    color: 'rgba(156, 165, 196, 0.5)',

    line: {

      color: 'rgba(156, 165, 196, 1)',

      width: 1,

    },

    symbol: 'circle',

    size: 20

  },

  hoverlabel: {

    bgcolor: 'black',

  }

};



var traceB = {

  type: "scatter",

  mode: "markers",

  x: [2011, 2012, 2013, 2014, 2015, 2016],

  y: [769, 755, 747, 770, 771, 781],

  name: 'Second Highest Marks',

  marker: {

    color: 'rgba(165, 196, 50, 0.5)',

    line: {

      color: 'rgba(165, 196, 50, 1)',

      width: 1,

    },

    symbol: 'circle',

    size: 20

  },

  hoverlabel: {

    bgcolor: 'black',

  }

};



var data = [traceA, traceB];



var layout = {

  title: 'Marks Obtained by Top Two Students',

  xaxis: {

    showgrid: false,

    showline: true,

    linecolor: 'rgb(200, 0, 0)',

    ticks: 'outside',

    tickcolor: 'rgb(200, 0, 0)',

    tickwidth: 4

  },

  legend: {

    bgcolor: 'white',

    borderwidth: 1,

    bordercolor: 'black',

    orientation: 'h',

    xanchor: 'center',

    x: 0.5,

    font: {

      size: 12,

    }

  },

  paper_bgcolor: 'rgb(255, 230, 255)',

  plot_bgcolor: 'rgb(255, 230, 255)'

};



Plotly.plot(dotDiv, data, layout);