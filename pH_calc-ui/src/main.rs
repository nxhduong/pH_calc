#![allow(non_snake_case, reason = "pH_min, pKa etc.")]

use cursive::{
    align::HAlign, 
    views::{
        Button, 
        Dialog, 
        LinearLayout, 
        TextContent, 
        TextView
    }, 
    Cursive, 
    CursiveExt
};
use cursive_table_view::{
    TableView,
    TableViewItem
};
use pH_calc_core::{
    types::AcidBase,
    calculator::compute_pH
};

#[derive(Copy, Clone, PartialEq, Eq, Hash)]
enum SpeciesColumn
{
    Name,
    Type,
    Conc,
    Consts
}

impl TableViewItem<SpeciesColumn> for AcidBase {

    fn to_column(&self, column: SpeciesColumn) -> String {
        match column {
            SpeciesColumn::Name => self.name.clone().unwrap_or("Unnamed".to_string()),
            SpeciesColumn::Type => String::from(if self.is_acidic() { "Acid" } else { "Base" }),
            SpeciesColumn::Conc => format!("{}", self.conc()),
            SpeciesColumn::Consts => format!("{:?}", self.dissoc_consts_acid())
        }
    }

    fn cmp(&self, other: &Self, column: SpeciesColumn) -> std::cmp::Ordering where Self: Sized {
        /* not needed 
        match column {
            SpeciesColumn::Name => self.name.cmp(&other.name),
            SpeciesColumn::Type => self.is_acidic().cmp(&other.is_acidic()),
            SpeciesColumn::Conc => self.conc().cmp(&other.conc())
            SpeciesColumn::Consts => self.dissoc_consts_acid().into_iter().cmp(&other.dissoc_consts_acid())
        }*/
        std::cmp::Ordering::Equal
    }
}

fn main()
{
    let mut species: Vec<AcidBase> = Vec::new();
    let mut root = Cursive::default();
    let mut result_text_label = TextView::new_with_content(TextContent::new("pH: x"));

    root.add_layer(
    Dialog::around(
        LinearLayout::horizontal()
            .child(
            TableView::<AcidBase, SpeciesColumn>::new()
                .column(SpeciesColumn::Name, "Name", |item| item)
                .column(SpeciesColumn::Type, "Acid/Base", |item| item)
                .column(SpeciesColumn::Conc, "Concentration", |item| item.align(HAlign::Center))
                .column(SpeciesColumn::Consts, "pKa Values", |item| item.align(HAlign::Right))
            )
            .child(
                LinearLayout::vertical()
                    .child(Button::new("Add", add_species))
                    .child(Button::new("Delete a species", del_species))
                    .child(Button::new("Calculate pH", calc_pH))
                    .child(result_text_label)
            )
        )
    );
    
    root.run();
}

fn add_species(cur: &mut Cursive)
{
    todo!()
}

fn del_species(cor: &mut Cursive)
{
    // By name
}

fn calc_pH(cur: &mut Cursive)
{
    todo!()
}