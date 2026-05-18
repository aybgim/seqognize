#[macro_use]
extern crate yew;

use yew::prelude::*;
use yew::services::{DialogService};
use seqognize::nt_aligner::{GlobalNtAligner, NtAlignmentConfig};
use seqognize::aligner::Aligner;
use std::num::ParseIntError;
use seqognize::element::Score;

struct Model {
    reference: String,
    subject: String,
    match_score: String,
    mismatch_score: String,
    alignment: String,
    score: String,
    parser: Parser,
}

enum Msg {
    SetReference(String),
    SetSubject(String),
    SetMatchScore(String),
    SetMismatchScore(String),
    Align,
}

impl Component for Model {
    type Message = Msg;
    type Properties = ();

    fn create(_: Self::Properties, _: ComponentLink<Self>) -> Self {
        Self {
            reference: "".to_string(),
            subject: "".to_string(),
            match_score: "1".to_string(),
            mismatch_score: "-1".to_string(),
            alignment: "".to_string(),
            score: "".to_string(),
            parser: Parser::new(),
        }
    }

    fn update(&mut self, msg: Self::Message) -> ShouldRender {
        match msg {
            Msg::Align => {
                let result = self.align().unwrap_or(AlignmentResult::empty());
                self.alignment = result.alignment;
                self.score = match result.score {
                    Some(s) => s.to_string(),
                    _ => "".to_string()
                };
            }
            Msg::SetSubject(value) => {
                self.reference = value;
                self.reset();
            }
            Msg::SetReference(value) => {
                self.subject = value;
                self.reset()
            }
            Msg::SetMatchScore(value) => {
                self.match_score = value;
                self.reset();
            }
            Msg::SetMismatchScore(value) => {
                self.mismatch_score = value;
                self.reset();
            }
        }
        true
    }
}

impl Model {
    fn reset(&mut self) {
        self.alignment = "".to_string();
        self.score = "".to_string();
    }
}

impl Renderable<Model> for Model {
    fn view(&self) -> Html<Self> {
        html! {
            <table>
                <tr>
                    <td>{"Reference:"}</td>
                    <td>
                        <input size="60", oninput=|e| Msg::SetReference(e.value),/>
                    </td>
                </tr>
                <tr>
                    <td>{"Subject:"}</td>
                    <td>
                        <input size="60", oninput=|e| Msg::SetSubject(e.value),/>
                    </td>
                </tr>
                <tr>
                    <td>{"Match score:"}</td>
                    <td>
                        <input type="numerical", size="5",
                            value={&self.match_score},
                            oninput=|e| Msg::SetMatchScore(e.value),
                        />
                    </td>
                </tr>
                <tr>
                    <td>{"Mismatch score:"}</td>
                    <td>
                        <input type="numerical", size="5",
                            value={&self.mismatch_score},
                            oninput=|e| Msg::SetMismatchScore(e.value),
                        />
                    </td>
                </tr>
                <tr>
                    <td align="left", >
                        <button onclick=|_| Msg::Align,>
                            {"Align"}
                        </button>
                    </td>
                    <td>
                        <textarea readonly="true", rows="3", cols="60",>
                            {&self.alignment}
                        </textarea>
                    </td>
                </tr>
                <tr>
                    <td>{"Alignment score:"}</td>
                    <td align="left", >
                        <input type="numerical", size="5", readonly="true",
                            value={&self.score},
                        />
                    </td>
                </tr>
            </table>
        }
    }
}

struct Parser {
    dialog: DialogService
}

impl Parser {
    fn new() -> Self {
        Parser { dialog: DialogService::new() }
    }

    fn parse(&mut self, value: &str) -> Result<Score, ParseIntError> {
        match value.parse::<Score>() {
            Ok(number) => Ok(number),
            Err(e) => {
                let msg = format!("Invalid number: {}", value);
                self.dialog.alert(&msg);
                return Err(e);
            }
        }
    }
}

fn main() {
    yew::start_app::<Model>();
}

impl Model {
    fn config(&mut self) -> Result<NtAlignmentConfig, ParseIntError> {
        Ok(NtAlignmentConfig::new(
            self.parser.parse(&self.match_score)?,
            self.parser.parse(&self.mismatch_score)?,
            -1i16,
            -1i16,
        ))
    }

    fn align(&mut self) -> Result<AlignmentResult, ParseIntError> {
        let config = self.config()?;
        let aligner = GlobalNtAligner { config };
        let alignment = match aligner.align(
            &self.subject.as_bytes(),
            &self.reference.as_bytes(),
        ) {
            Ok(a) => a,
            Err(_) => return Ok(AlignmentResult::empty()),
        };
        let aligned_sequences = alignment.aligned_sequences();
        let alignment_str = format!("{}\n{}\n{}", aligned_sequences.0, aligned_sequences.1, aligned_sequences.2);
        Ok(AlignmentResult::of(alignment_str, alignment.score))
    }
}

struct AlignmentResult {
    alignment: String,
    score: Option<Score>,
}

impl AlignmentResult {
    fn of(alignment: String, score: Score) -> Self {
        AlignmentResult { alignment, score: Some(score) }
    }

    fn empty() -> Self {
        AlignmentResult { alignment: "".to_string(), score: None }
    }
}

