#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
#########################################################################
# Author: bqxiao
# File: draw_html.py
# Created time: 2025/06/19 16:10:14
# E-mail: xiao.benqi@qq.com
#########################################################################
# Usage: python draw_html.py -h


import matplotlib.pyplot as plt
import json
import argparse
import networkx as nx

helptext=r"""
例子：
python draw_html.py -i A01.num -a A01.txt -o A01.html

输入的文件格式为(这两个文件每个元素必须一一对应):

--input(必须) :
1	1	1	1	1	934	1	1	1	1	1	1	1	1
2	2	2	2	2	933	2	2	2	2	2	2	2	2
3	3	3	3	3	932	3	3	3	3	3	3	3	3
4	4	4	4	4	-	4	4	4	4	4	4	4	4
5	5	5	5	5	-	5	5	5	5	5	5	5	5
6	6	6	6	6	-	6	6	6	6	6	6	6	6

--annotation(非必须) 若是显示多行信息,请使用 ";" 作为换行符:
spe1;1;gene1;range:0-545455	spe2;1	spe3;1	spe4;1	spe5;1	spe6;934	spe7;1	spe8;1	spe9;1	spe10;1	spe11;1	spe12;1	spe13;1	spe14;1
spe1;2;gene2;range:545455-545700	spe2;2	spe3;2	spe4;2	spe5;2	spe6;933	spe7;2	spe8;2	spe9;2	spe10;2	spe11;2	spe12;2	spe13;2	spe14;2
spe1;3;gene3;range:545600-545900	spe2;3	spe3;3	spe4;3	spe5;3	spe6;932	spe7;3	spe8;3	spe9;3	spe10;3	spe11;3	spe12;3	spe13;3	spe14;3
spe1;4;	spe2;4	spe3;4	spe4;4	spe5;4	spe6;-	spe7;4	spe8;4	spe9;4	spe10;4	spe11;4	spe12;4	spe13;4	spe14;4
spe1;5;	spe2;5	spe3;5	spe4;5	spe5;5	spe6;-	spe7;5	spe8;5	spe9;5	spe10;5	spe11;5	spe12;5	spe13;5	spe14;5
spe1;6	spe2;6	spe3;6	spe4;6	spe5;6	spe6;-	spe7;6	spe8;6	spe9;6	spe10;6	spe11;6	spe12;6	spe13;6	spe14;6
"""

parser=argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-i","--input",type=str,required=True,help="geneindex num file(注释为基因顺序)")
parser.add_argument("-a","--annotation",type=str,required=True,help="geneindex txt file(注释文件,必须与num文件一一对应)")
parser.add_argument("-o","--html",type=str,required=True,help="output html file")
parser.add_argument("-l","--lazy",action="store_false",help="是否使用惰性加载,默认为 True")
parser.add_argument("-r","--range",help="显示范围,格式:起始-结束,例如:100-200,默认为全部显示")
parser.add_argument("-s","--species",type=str,help="指定物种,格式:1-15,17,18 (表示显示第0到15个物种,以及第17和18个物种,逗号分隔物种),默认为全部物种")
args=parser.parse_args()

############################ 网络图基础配置信息 ############################
YGAP=40 # 元素之间的间距 y(竖向间隔)
XGAP=70 # 元素之间的间距 x(横向间隔)

noAdjacentCurve=40 # 非相邻元素之间的曲线弧度
noAdjacentCurveColor="#888" # 非相邻元素之间的曲线颜色
noAdjacentCurveWidth=2 # 非相邻元素之间的曲线宽度
AdjacentCurve=0 # 相邻元素之间的曲线弧度

CurveWidth=2 # 每个物种基因邻接网络图：曲线宽度
ancestorCurveWidthMax=5 # pan基因邻接网络图：曲线宽度最大值,这里我们设想边的权值与宽度正相关


CurveColor="#888" # 曲线颜色
ancestorCurveColor="#888" # 祖先曲线颜色

InversionCurveColor="#FF0000" # 倒位曲线颜色
InversionCurveWidth=5 # 倒位曲线宽度

#nodeShape="circle" # 默认节点形状 circle / roundrectangle
nodewidth=40 # 若是使用矩形节点，设置节点宽度
nodeheight =40 # 若是使用矩形节点，设置节点高度

targetArrowCurveShape = "vee" # 箭头形状 : triangle, vee, circle, diamond
############################ 网络图基础配置信息 ############################


############################ 函数功能 ############################
def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % tuple(int(255*x) for x in rgb)
tab20_rgb = plt.get_cmap('tab20').colors
css_colors = [rgb_to_hex(c) for c in tab20_rgb]

def parse_range(range_str,max_length:int)-> tuple[int,int] :
    """解析范围字符串，返回起始和结束索引的元组"""
    try:
        start, end = map(int, range_str.split('-'))
        if start > end :
            start, end = end, start
        start=start-1  # 转换为0基索引
        if start < 0:
            raise ValueError
        if end > max_length:
            end = max_length
        return start , end  # 转换为0基索引
    except Exception:
        raise argparse.ArgumentTypeError("范围格式错误，应为 '起始-结束'，例如 '100-200'，且起始应小于等于结束")

def parse_species(species_str,max_species:int)-> list[int] :
    """解析物种字符串，返回物种索引的列表"""
    species_set = set()
    spe_l=species_str.split(',')
    for part in spe_l:
        if "-" in part:
            try:
                start, end = map(int, part.split('-'))
                if start > end :
                    start, end = end, start
                if start < 1 or end > max_species:
                    raise ValueError(f"物种索引应在[{start},{end}] [1,{max_species}] 之间")
                start=start-1  # 转换为0基索引
                end=end      # 转换为0基索引
                species_set.update(range(start, end))
            except Exception:
                raise argparse.ArgumentTypeError("物种格式错误，应为 '起始-结束'，例如 '1-15' ")
        else:
            try:
                idx = int(part)
                if idx < 1 or idx > max_species:
                    raise ValueError(f"物种索引应在 [1,{max_species}] 之间")
                idx=idx-1  # 转换为0基索引
                species_set.add(idx)
            except Exception:
                raise argparse.ArgumentTypeError("物种格式错误，应为整数或 '起始-结束'，例如 '1,2,3'")
    spe_l=list(species_set)
    spe_l.sort()
    return spe_l

############################ 函数功能 ############################

infile=args.input
outfile=args.html
lazy=args.lazy
annotationfile=args.annotation
species_list=args.species


data=[]
for i in open(infile):
    i=i.strip().split("\t")
    line = [int(z) if z != "-" else -9 for z in i]
    data.append(line)
data_len=len(data) if data else 0
data_species_num=len(data[0]) if data else 1
data=list(zip(*data))

annotation_data=[]
for i in open(annotationfile):
    i=i.strip().split("\t")
    annotation_data.append(i)
annotation_data_len=len(annotation_data) if annotation_data else 0
annotation_data_species_num=len(annotation_data[0]) if annotation_data else 1
annotation_data=list(zip(*annotation_data))

# 检查注释文件和基因位置文件是否匹配，目的就是使得 .num文件 的ID和 注释文件 的ID可以对应上
if annotation_data_len != data_len:
    raise ValueError(f"注释文件行数 {annotation_data_len} 和 基因位置行数 {data_len} 不匹配，请检查输入文件")
if annotation_data_species_num != data_species_num:
    raise ValueError(f"注释文件列数 {annotation_data_species_num} 和 基因位置列数 {data_species_num} 不匹配，程序继续运行，但请检查输入文件")



if args.range:
    min_start,max_end=parse_range(args.range,data_len)
else:
    min_start=0
    max_end=data_len
if args.species:
    species_idx=parse_species(args.species,data_species_num)
    data=[data[i] for i in species_idx]
    annotation_data=[annotation_data[i] for i in species_idx]
    data_species_num=len(data)
data=[i[min_start:max_end] for i in data]
annotation_data=[i[min_start:max_end] for i in annotation_data]

# 节点ID到注释的映射，这里将有逗号替换为换行符
id2annotation={}
for i_index,i in enumerate(annotation_data):
    for j_index,j in enumerate(i):
        id=f"{i_index}_{j_index}"
        j=j.replace(",","<br>")
        id2annotation[id]=j
id2num={}
for i_index,i in enumerate(data):
    for j_index,j in enumerate(i):
        id=f"{i_index}_{j_index}"
        id2num[id]=j

G=nx.DiGraph()
nodes_edges=[]

# 给每一个物种添加节点
y_position=YGAP*2
for i_index,i in enumerate(annotation_data):
    y_position+=YGAP
    x_position=0
    for j_index,j in enumerate(i):
        x_position+=XGAP
        id=f"{i_index}_{j_index}"
        color=css_colors[j_index%20]
        
        # 添加节点
        this_annotation=id2annotation.get(id,id)
        this_gene_num=this_annotation.count("<br>")+1
        this_lable=id2num.get(id,j)
        if this_lable==-9 :
            if this_annotation=="-":
                continue
            else:
                this_lable=""
        if this_gene_num==1:
            thisNodeClass="single"
        else:
            thisNodeClass="multiple"
        nodes_edges.append({"data":{"id":id,
                                    "color":color,
                                    "label":this_lable,
                                    "showtxt":this_annotation},
                            "position":{"x":x_position,
                                        "y":y_position},
                            "classes":thisNodeClass
                            })

# 给每一个物种添加边
for i_index,i in enumerate(data):
    # 记录每个节点的index，避免混淆节点id
    num2index={}
    for j_index,j in enumerate(i):
        if j==-9:
            continue
        num2index[j]=j_index
    
    noNoneList=[j for j in i if j != -9]
    
    # 判断两个节点在祖先中是否相邻，若是相邻则添加直线边，否则就是弧线边
    num2NoGapindex={}
    for j_index,j in enumerate(noNoneList):
        num2NoGapindex[j]=j_index

    # 曲线和直线边
    noNoneList.sort()
    for z in range(len(noNoneList)-1):
        thisnode=noNoneList[z]
        index_thisnode=num2index[noNoneList[z]]
        nextnode=noNoneList[z+1]
        index_nextnode=num2index[noNoneList[z+1]]
        if G.has_edge(index_thisnode,index_nextnode):
            G[index_thisnode][index_nextnode]["weight"]+=1
        else:
            G.add_edge(index_thisnode,index_nextnode,weight=1)
        id=f"edge{i_index}_{index_thisnode}_{index_nextnode}"
        source=f"{i_index}_{index_thisnode}"
        target=f"{i_index}_{index_nextnode}"
        color=CurveColor
        width=CurveWidth
        # 判断两个节点在祖先中书否相邻
        if abs(num2NoGapindex[thisnode]-num2NoGapindex[nextnode])==1:
            curve=AdjacentCurve
            if num2NoGapindex[thisnode]<num2NoGapindex[nextnode]:
                pass
            else:
                color=InversionCurveColor
                width=InversionCurveWidth
        else:
            curve=noAdjacentCurve
            color=noAdjacentCurveColor
            width=noAdjacentCurveWidth
        label=""
        # 添加边
        nodes_edges.append({"data":{"id":id,
                                    "source":source,
                                    "width":width,
                                    "target":target,
                                    "curve":curve,
                                    "label":label,
                                    "color":color}})

# pan 基因集节点
y_position=0
pan_genes=list(range(min_start,max_end+1))
x_position=0
for kz,SGid in enumerate(pan_genes):
    x_position+=XGAP
    id=f"anc_{SGid}"
    color=css_colors[kz%20]
    label=SGid
    nodes_edges.append({"data":{"id":id,
                                "color":color,
                                "label":label,
                                "showtxt":f"SG{SGid:07d}"},
                        "position":{"x":x_position,
                                    "y":y_position},
                        "classes":"single"
                                    })
# pan 基因集边
for edge in G.edges(data=True):
    thisnode=edge[0]+min_start
    nextnode=edge[1]+min_start
    weight=edge[2]["weight"]
    id=f"edge_anc_{thisnode}_{nextnode}"
    source=f"anc_{thisnode}"
    target=f"anc_{nextnode}"
    color=ancestorCurveColor
    width=weight/data_species_num*ancestorCurveWidthMax
    curve=noAdjacentCurve
    #text=weight
    text=""
    nodes_edges.append({"data":{"id":id,
                                "source":source,
                                "width":width,
                                "target":target,
                                "curve":curve,
                                "label":text}})


nodes_edges=[str(i) for i in nodes_edges]
nodeEdges="["+',\n'.join(nodes_edges)+"]"

# CSS 样式
style = """

/* CSS 样式 */
html, body {
    margin: 0;
    padding: 0;
    width: 100%;
    height: 100%;
    overflow: hidden;
}

/* Cytoscape 布局样式 */
#cy {
    width: 100%;
    height: 100%;
    display: block;
}

/* tooltip 样式 */
.cy-tooltip {
  position: absolute;
  background: rgba(0,0,0,0.75);
  color: white;
  padding: 6px 10px;
  border-radius: 4px;
  font-size: 13px;
  pointer-events: none;
  z-index: 9999;
  display: none;
  white-space: pre-line; /* 支持换行 */
}

/* 复制成功提示样式 */
#copyStatus {
  position: fixed;       /* 固定定位，脱离文档流，始终相对于视口 */
  top: 10px; right: 10px; /* 定位在右上角（距顶部和右侧各10px） */
  background: #498db4;   /* 背景色：浅蓝色 */
  color: white;          /* 文字颜色：白色 */
  padding: 8px 12px;     /* 内边距：上下8px，左右12px */
  border-radius: 4px;    /* 圆角：4px */
  font-size: 14px;       /* 字体大小：14px */
  opacity: 0;            /* 默认完全透明（隐藏状态） */
  transition: opacity 0.3s ease; /* 透明度变化过渡动画：0.3秒缓动 */
  pointer-events: none;  /* 默认不响应鼠标事件（避免遮挡下方元素） */
  z-index: 10000;        /* 层级极高，确保悬浮在其他元素上方 */
}

/* 显示复制成功提示 */
#copyStatus.show {
  opacity: 1;            /* 完全不透明（显示状态） */
  pointer-events: auto;  /* 允许响应鼠标事件（如需交互） */
}
"""

# Cytoscape 样式
cy_style = [
    # 节点样式
    {
        "selector": "node.single",
        "style": {
            "background-color": "data(color)",
            "label": "data(label)",
            "color": "#fff",#白色字体
            "text-valign": "center",
            "text-halign": "center",
            "font-size": "12px",
            "font-family": "Arial",
        }
    },
    # 节点样式
    {
        "selector": "node.multiple",
        "style": {
            "background-color": "data(color)",
            "shape": "roundrectangle",
            "width": nodewidth,
            "height": nodeheight,
            "label": "data(label)",
            "color": "#fff",#白色字体
            "text-valign": "center",
            "text-halign": "center",
            "font-size": "12px",
            "font-family": "Arial",
        }
    },
    # 边样式
    {
        "selector": "edge",
        "style": {
            "line-color": "data(color)",
            'target-arrow-color': 'data(color)',
            "label": "data(label)",
            "width": "data(width)",
            "target-arrow-shape": targetArrowCurveShape ,#箭头样式
            "curve-style": "bezier" 
        }
    }
]

tooltip=r"""
const tooltip = document.getElementById('tooltip');
const copyStatus = document.getElementById('copyStatus');

// 显示 tooltip
cy.on('mouseover', 'node', (event) => {
  const node = event.target;
  const description = node.data('showtxt') || '';
  tooltip.innerHTML = description;
  tooltip.style.display = 'block';
});

// 隐藏 tooltip
cy.on('mouseout', 'node', () => {tooltip.style.display = 'none';});

// tooltip 跟随鼠标
cy.on('mousemove', (event) => {
  tooltip.style.left = event.originalEvent.pageX + 10 + 'px';
  tooltip.style.top = event.originalEvent.pageY + 10 + 'px';
});

// 点击节点复制 tooltip 内容到剪贴板
cy.on('tap', 'node', async (event) => {
  const node = event.target;
  let text = node.data('showtxt') || '';

  // 去掉 HTML 标签，只保留纯文本（否则复制时会带上 <br> 标签）
  const newStr = text.replace(/<br>/g, "\n");
  try {
    // 使用 Clipboard API 复制
    await navigator.clipboard.writeText(newStr);

    // 复制成功提示词
    copyStatus.textContent = newStr + ' 复制成功 ';

    // 显示复制成功提示
    copyStatus.classList.add('show');
    setTimeout(() => {copyStatus.classList.remove('show');}, 1000);
  } catch (err) {
    alert('复制失败');
  }
});
"""

# 延迟加载脚本
LazyLoadingScript = """
// 获取当前视窗中可见的元素
    function getVisibleElements() {
        const pan = cy.pan();
        const zoom = cy.zoom();

        const viewport = {
            left: -pan.x / zoom,
            top: -pan.y / zoom,
            right: (-pan.x + window.innerWidth) / zoom,
            bottom: (-pan.y + window.innerHeight) / zoom
        };

        const visibleNodes = elements.filter(el => {
            return el.position &&
                el.position.x >= viewport.left &&
                el.position.x <= viewport.right &&
                el.position.y >= viewport.top &&
                el.position.y <= viewport.bottom;
        });

        const visibleIds = new Set(visibleNodes.map(n => n.data.id));

        const visibleEdges = elements.filter(el => {
            return el.data && el.data.source && el.data.target &&
                visibleIds.has(el.data.source) &&
                visibleIds.has(el.data.target);
        });

        return [...visibleNodes, ...visibleEdges];
    }

    function updateVisibleElements() {
        const visible = getVisibleElements();
        cy.elements().remove();
        cy.add(visible);
        cy.layout({ name: 'preset' }).run();

        // 设置曲线
        cy.edges().forEach(edge => {
            const curve = edge.data('curve');
            if (curve === 0) {
                edge.style({ 'curve-style': 'straight' });
            } else {
                edge.style({
                    'curve-style': 'unbundled-bezier',
                    'control-point-distances': [curve],
                    'control-point-weights': [0.5]
                });
            }
        });
    }

    // 初次渲染
    updateVisibleElements();

    // 添加拖动/缩放监听，节流处理更新
    let updateTimer = null;
    cy.on('pan zoom', () => {
        if (updateTimer) return;
        updateTimer = setTimeout(() => {
            updateVisibleElements();
            updateTimer = null;
        }, 200);
    });
"""

NOTLazyLoadingScript=f"""
// 设置曲线
cy.edges().forEach(edge => {{
    const curve = edge.data('curve');
    if (curve === {AdjacentCurve}) {{
        edge.style({{ 'curve-style': 'straight' }});
    }} else {{
        edge.style({{
            'curve-style': 'unbundled-bezier',
            'control-point-distances': [curve],
            'control-point-weights': [0.5]
        }});
    }}
}});
"""

if not lazy:
    # 不开启延迟加载
    # 使用 --lazy
    LazyLoadingScript=NOTLazyLoadingScript
    elementsTextInCytoscape = "elements: elements"
    print("The web page does not use lazy loading")
else:
    # 开启延迟加载
    # 不使用 --lazy
    elementsTextInCytoscape = ""
    print("Lazy loading of web pages")

# 输出 HTML
html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset='utf-8'>
    <title></title>
    <script src='https://unpkg.com/cytoscape@3.26.0/dist/cytoscape.min.js'></script>
    <style>
        {style}
    </style>
</head>
<body>
    <div id='cy'></div>
    <!-- tooltip -->
    <div id="tooltip" class="cy-tooltip"></div>
    <!-- 复制提示 -->
    <div id="copyStatus"></div>
    <script>
    const elements = {nodeEdges}
        const cy = cytoscape({{
            container: document.getElementById('cy'),
            layout: {{ name: 'preset' }},
            style: {json.dumps(cy_style)},
            {elementsTextInCytoscape}
            
        }});

        {LazyLoadingScript}
        {tooltip}
    </script>
</body>
</html>
"""

# 写入 HTML 文件
with open(outfile, "w", encoding="utf-8") as f:
    f.write(html)

print(f"HTML 文件已生成: {outfile}")
