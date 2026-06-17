const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  Header, Footer, AlignmentType, LevelFormat, HeadingLevel, BorderStyle,
  WidthType, ShadingType, PageNumber, PageBreak, TabStopType, TabStopPosition,
} = require("docx");

// ---- Load data: real pipeline output if available, else 10 placeholders ----
const DATA_PATH = process.argv[2] || path.join(__dirname, "calibration_run.json");

let entries;
if (fs.existsSync(DATA_PATH)) {
  const raw = JSON.parse(fs.readFileSync(DATA_PATH, "utf8"));
  entries = raw.results.map((r, i) => ({
    index: r.index || i + 1,
    topic: r.topic || "",
    hypothesis: r.hypothesis?.hypothesis || "[hypothesis text unavailable]",
    experimental_approach: r.hypothesis?.experimental_approach || "",
    fly_lines: (r.hypothesis?.fly_lines_needed || [])
      .map(fl => `${fl.line || "?"} — ${fl.purpose || ""} (${fl.source || "source unspecified"})`),
    expected_outcomes: r.hypothesis?.expected_outcomes || "",
    novelty_claim: r.hypothesis?.novelty_claim || "",
  }));
} else {
  console.log(`No data file found at ${DATA_PATH} — generating a 10-slot blank template instead.`);
  entries = Array.from({ length: 10 }, (_, i) => ({
    index: i + 1,
    topic: "[research topic — to be filled in from pipeline output]",
    hypothesis: "[hypothesis text — to be filled in from pipeline output]",
    experimental_approach: "[experimental approach — to be filled in]",
    fly_lines: [],
    expected_outcomes: "",
    novelty_claim: "",
  }));
}

// ---- Style helpers ----
const FONT = "Arial";
const ACCENT = "2E5C8A";
const LIGHT_FILL = "EAF1F8";
const border = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: border, bottom: border, left: border, right: border };

const PAGE_WIDTH = 12240;
const PAGE_HEIGHT = 15840;
const MARGIN = 1080; // 0.75"
const CONTENT_WIDTH = PAGE_WIDTH - MARGIN * 2;

function heading(text, level = HeadingLevel.HEADING_1) {
  return new Paragraph({ heading: level, children: [new TextRun(text)] });
}

function bodyText(text, opts = {}) {
  return new Paragraph({
    spacing: { after: 120 },
    children: [new TextRun({ text, ...opts })],
  });
}

function labeledBlock(label, text) {
  return new Paragraph({
    spacing: { after: 100 },
    children: [
      new TextRun({ text: `${label}: `, bold: true }),
      new TextRun({ text: text || "—" }),
    ],
  });
}

// 1-5 score row: label cell + 5 blank circle cells + a free "score:" line cell
function scoreRow(dimensionLabel, dimensionDesc) {
  const labelCellWidth = 3400;
  const scaleCellWidth = Math.floor((CONTENT_WIDTH - labelCellWidth) / 5);
  const lastCellWidth = CONTENT_WIDTH - labelCellWidth - scaleCellWidth * 4;

  const scaleCells = [1, 2, 3, 4, 5].map((n, i) =>
    new TableCell({
      borders,
      width: { size: i === 4 ? lastCellWidth : scaleCellWidth, type: WidthType.DXA },
      verticalAlign: "center",
      margins: { top: 100, bottom: 100, left: 60, right: 60 },
      children: [
        new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [new TextRun({ text: "◯", size: 36, color: "999999" })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { before: 20 },
          children: [new TextRun({ text: String(n), size: 18, color: "777777" })],
        }),
      ],
    })
  );

  return new TableRow({
    children: [
      new TableCell({
        borders,
        width: { size: labelCellWidth, type: WidthType.DXA },
        shading: { fill: LIGHT_FILL, type: ShadingType.CLEAR },
        margins: { top: 100, bottom: 100, left: 120, right: 120 },
        children: [
          new Paragraph({ children: [new TextRun({ text: dimensionLabel, bold: true, size: 21 })] }),
          new Paragraph({ children: [new TextRun({ text: dimensionDesc, size: 17, italics: true, color: "555555" })] }),
        ],
      }),
      ...scaleCells,
    ],
  });
}

function scaleHeaderRow() {
  const labelCellWidth = 3400;
  const scaleCellWidth = Math.floor((CONTENT_WIDTH - labelCellWidth) / 5);
  const lastCellWidth = CONTENT_WIDTH - labelCellWidth - scaleCellWidth * 4;

  const headers = ["1\nPoor", "2", "3\nAdequate", "4", "5\nExcellent"];
  const cells = headers.map((h, i) =>
    new TableCell({
      borders,
      width: { size: i === 4 ? lastCellWidth : scaleCellWidth, type: WidthType.DXA },
      shading: { fill: ACCENT, type: ShadingType.CLEAR },
      margins: { top: 80, bottom: 80, left: 40, right: 40 },
      children: h.split("\n").map(line =>
        new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [new TextRun({ text: line, size: 16, bold: true, color: "FFFFFF" })],
        })
      ),
    })
  );

  return new TableRow({
    children: [
      new TableCell({
        borders,
        width: { size: labelCellWidth, type: WidthType.DXA },
        shading: { fill: ACCENT, type: ShadingType.CLEAR },
        margins: { top: 80, bottom: 80, left: 120, right: 120 },
        children: [new Paragraph({ children: [new TextRun({ text: "Dimension", bold: true, size: 18, color: "FFFFFF" })] })],
      }),
      ...cells,
    ],
  });
}

function notesRow() {
  return new TableRow({
    children: [
      new TableCell({
        borders,
        columnSpan: 6,
        width: { size: CONTENT_WIDTH, type: WidthType.DXA },
        margins: { top: 120, bottom: 200, left: 120, right: 120 },
        children: [
          new Paragraph({ children: [new TextRun({ text: "Notes / comments:", bold: true, size: 19 })] }),
          new Paragraph({ children: [new TextRun({ text: " ", size: 19 })], border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "AAAAAA" } }, spacing: { before: 200, after: 200 } }),
          new Paragraph({ children: [new TextRun({ text: " ", size: 19 })], border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "AAAAAA" } }, spacing: { after: 200 } }),
          new Paragraph({ children: [new TextRun({ text: " ", size: 19 })], border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "AAAAAA" } }, spacing: { after: 100 } }),
        ],
      }),
    ],
  });
}

function scoreTable() {
  return new Table({
    width: { size: CONTENT_WIDTH, type: WidthType.DXA },
    rows: [
      scaleHeaderRow(),
      scoreRow("Hallucination Defense", "Avoids fabricated genes/citations/facts; correctly separates established fact from speculation"),
      scoreRow("Genetic Feasibility", "Fly lines, crosses, and tools are real and executable on a realistic timeline"),
      scoreRow("Control Reasoning", "Appropriate genetic/statistical/environmental controls are specified and justified"),
      notesRow(),
    ],
  });
}

function hypothesisBlock(entry) {
  const flyLinesText = entry.fly_lines.length
    ? entry.fly_lines.join("\n")
    : "—";

  return [
    new Paragraph({
      heading: HeadingLevel.HEADING_2,
      spacing: { before: 360, after: 120 },
      children: [new TextRun(`Hypothesis ${entry.index}`)],
    }),
    labeledBlock("Topic", entry.topic),
    labeledBlock("Hypothesis", entry.hypothesis),
    labeledBlock("Experimental approach", entry.experimental_approach),
    labeledBlock("Fly lines needed", flyLinesText),
    labeledBlock("Expected outcomes", entry.expected_outcomes),
    labeledBlock("Stated novelty / gap claim", entry.novelty_claim),
    new Paragraph({ spacing: { before: 160, after: 80 }, children: [new TextRun({ text: "Grader scores (1–5 scale):", bold: true, size: 21 })] }),
    scoreTable(),
  ];
}

// ---- Cover / instructions page ----
const coverChildren = [
  new Paragraph({
    heading: HeadingLevel.TITLE,
    spacing: { after: 200 },
    children: [new TextRun("Drosophila Hypothesis Pipeline — Calibration Scoring Sheet")],
  }),
  bodyText(
    "This sheet contains hypotheses generated by an AI research-assistance pipeline for Drosophila melanogaster cancer and developmental biology research. Each hypothesis was produced from a real literature synthesis and includes a proposed experimental design. Please score each hypothesis independently across the three dimensions below using a 1–5 scale. There is no pass/fail threshold — please rate quality as you see it."
  ),
  heading("Scoring dimensions", HeadingLevel.HEADING_2),
  new Paragraph({
    spacing: { after: 100 },
    children: [
      new TextRun({ text: "Hallucination Defense (1–5): ", bold: true }),
      new TextRun("Does the hypothesis avoid asserting fabricated genes, citations, or facts? Does it correctly distinguish established biology from speculation, rather than overclaiming certainty?"),
    ],
  }),
  new Paragraph({
    spacing: { after: 100 },
    children: [
      new TextRun({ text: "Genetic Feasibility (1–5): ", bold: true }),
      new TextRun("Are the proposed fly lines, crosses, and genetic tools real, obtainable, and executable by a grad student on a realistic timeline?"),
    ],
  }),
  new Paragraph({
    spacing: { after: 200 },
    children: [
      new TextRun({ text: "Control Reasoning (1–5): ", bold: true }),
      new TextRun("Are appropriate genetic, statistical, and environmental controls specified and justified, such that the result would be interpretable as described?"),
    ],
  }),
  heading("Scale reference", HeadingLevel.HEADING_2),
  bodyText("1 = Poor — significant, substantive problems on this dimension.", { size: 21 }),
  bodyText("2 = Below average — notable problems, but not severe.", { size: 21 }),
  bodyText("3 = Adequate — workable, with some room for improvement.", { size: 21 }),
  bodyText("4 = Good — minor issues only.", { size: 21 }),
  bodyText("5 = Excellent — no meaningful issues on this dimension.", { size: 21 }),
  new Paragraph({ children: [new PageBreak()] }),
];

// ---- Assemble document ----
const doc = new Document({
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 30, bold: true, font: FONT, color: ACCENT },
        paragraph: { spacing: { before: 240, after: 200 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 25, bold: true, font: FONT, color: ACCENT },
        paragraph: { spacing: { before: 280, after: 140 }, outlineLevel: 1 } },
      { id: "Title", name: "Title", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 36, bold: true, font: FONT, color: ACCENT },
        paragraph: { spacing: { after: 240 } } },
    ],
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: PAGE_WIDTH, height: PAGE_HEIGHT },
          margin: { top: MARGIN, right: MARGIN, bottom: MARGIN, left: MARGIN },
        },
      },
      headers: {
        default: new Header({
          children: [
            new Paragraph({
              alignment: AlignmentType.RIGHT,
              children: [new TextRun({ text: "Drosophila AI Pipeline — Calibration Round", size: 16, color: "888888" })],
            }),
          ],
        }),
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.CENTER,
              children: [
                new TextRun({ text: "Grader name: __________________________      Page ", size: 16, color: "888888" }),
                new TextRun({ children: [PageNumber.CURRENT], size: 16, color: "888888" }),
              ],
            }),
          ],
        }),
      },
      children: [
        ...coverChildren,
        ...entries.flatMap((e, i) => [
          ...hypothesisBlock(e),
          ...(i < entries.length - 1 ? [new Paragraph({ children: [new PageBreak()] })] : []),
        ]),
      ],
    },
  ],
});

Packer.toBuffer(doc).then(buffer => {
  const outPath = path.join(__dirname, "scoring_sheet.docx");
  fs.writeFileSync(outPath, buffer);
  console.log(`Wrote ${outPath}`);
});
