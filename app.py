from flask import Flask, render_template, request, jsonify
import pandas as pd
import plotly.express as px

app = Flask(__name__)

filtered_df = pd.read_csv('filtered_df.csv')

@app.route("/")
def home():
    return render_template("home.html")

@app.route("/index")
def index():
    primary_sites = filtered_df['PRIMARY_SITE'].unique()
    return render_template("index.html", primary_sites=primary_sites)

@app.route("/organ/<string:name>")
def get_organ_info(name):
    if name not in filtered_df['PRIMARY_SITE'].unique():
        return jsonify({"error": "Invalid organ name"}), 404

    organ_data = filtered_df[filtered_df['PRIMARY_SITE'] == name]
    cancer_type_counts = organ_data['CANCER_TYPE'].value_counts().to_dict()

    return jsonify({"organ": name, "cancer_type_counts": cancer_type_counts})

@app.route("/info", methods=["GET"])
def organ_info():
    organ_name = request.args.get("organ")

    if organ_name is None or organ_name not in filtered_df['PRIMARY_SITE'].unique():
        return render_template("error.html", error="Invalid organ name or no organ selected", link="/")

    organ_data = filtered_df[filtered_df['PRIMARY_SITE'] == organ_name]
    cancer_type_counts = organ_data['CANCER_TYPE'].value_counts().to_dict()

    return render_template("info.html", organ=organ_name, cancer_type_counts=cancer_type_counts, link="/")

@app.route("/substitution")
def substitution_matrix():
    # Matrix for base_allele and mutant_allele columns
    allele_matrix = pd.crosstab(filtered_df['BASE_ALLELE'], filtered_df['MUTANT_ALLELE'], margins=True, margins_name='Total')

    fig1 = px.imshow(allele_matrix, labels=dict(x="Mutant Allele", y="Wild Type Allele"),
                     x=allele_matrix.columns, y=allele_matrix.index,
                     title='Substitution Matrix (Allele)', color_continuous_scale='YlGnBu')

    fig1.update_layout(annotations=[
        dict(x=i, y=j, text=str(allele_matrix.iloc[j, i]),
             showarrow=False, font=dict(color='black'))
        for i in range(len(allele_matrix.columns))
        for j in range(len(allele_matrix.index))
    ])

    fig1.update_layout(xaxis=dict(tickvals=list(range(len(allele_matrix.columns))),
                                  ticktext=allele_matrix.columns),
                        yaxis=dict(tickvals=list(range(len(allele_matrix.index))),
                                  ticktext=allele_matrix.index))
    fig1.update_layout(width=800, height=400)

    # Amino Acid Matrix
    substitution_matrix = pd.crosstab(filtered_df['WT_AA_3'], filtered_df['MT_AA_3'], margins=True, margins_name='Total')

    fig2 = px.imshow(substitution_matrix, labels=dict(x="Mutant Amino Acid", y="Wild Type Amino Acid"),
                     x=substitution_matrix.columns, y=substitution_matrix.index,
                     title='Substitution Matrix (Amino Acid)', color_continuous_scale='YlGnBu')

    fig2.update_layout(annotations=[
        dict(x=i, y=j, text=str(substitution_matrix.iloc[j, i]),
             showarrow=False, font=dict(color='black'))
        for i in range(len(substitution_matrix.columns))
        for j in range(len(substitution_matrix.index))
    ])

    fig2.update_layout(xaxis=dict(tickvals=list(range(len(substitution_matrix.columns))),
                                  ticktext=substitution_matrix.columns),
                        yaxis=dict(tickvals=list(range(len(substitution_matrix.index))),
                                  ticktext=substitution_matrix.index))
    fig2.update_layout(width=1200, height=800)

    return render_template("substitution.html", plot1=fig1.to_html(full_html=False), plot2=fig2.to_html(full_html=False))

@app.route("/predictor", methods=["GET", "POST"])
def predictor():
    if request.method == "POST":
        base_allele = request.form.get("base_allele")
        mutant_allele = request.form.get("mutant_allele")
        primary_site = request.form.get("primary_site")

        if base_allele is not None and mutant_allele is not None and primary_site is not None:
            prediction_data = filtered_df[
                (filtered_df['BASE_ALLELE'] == base_allele) &
                (filtered_df['MUTANT_ALLELE'] == mutant_allele) &
                (filtered_df['PRIMARY_SITE'] == primary_site)
            ]

            cancer_type_probabilities = (
                prediction_data['CANCER_TYPE'].value_counts(normalize=True) * 100
            ).round(2).to_dict()

            return render_template("result.html", base_allele=base_allele, mutant_allele=mutant_allele,
                                   primary_site=primary_site, cancer_type_probabilities=cancer_type_probabilities)

    # If the form is not submitted or there's an error, render the predictor page with the form
    base_alleles = filtered_df['BASE_ALLELE'].unique()
    mutant_alleles = filtered_df['MUTANT_ALLELE'].unique()
    primary_sites = filtered_df['PRIMARY_SITE'].unique()

    return render_template("predictor.html", base_alleles=base_alleles, mutant_alleles=mutant_alleles, primary_sites=primary_sites)


@app.route("/histo")
def histo():
    tumour_types = filtered_df['TUMOUR_ORIGIN'].unique()

    figs = []
    for tumour_type in tumour_types:
        tumour_data = filtered_df[filtered_df['TUMOUR_ORIGIN'] == tumour_type]

        fig = px.bar(tumour_data, x='CANCER_TYPE', title=f'Cancer Types for {tumour_type}',
                     labels={'CANCER_TYPE': 'Cancer Type', 'count': 'Count'},
                     color='CANCER_TYPE')

        unique_cancer_types = tumour_data['CANCER_TYPE'].unique()
        for cancer_type in unique_cancer_types:
            count = tumour_data[tumour_data['CANCER_TYPE'] == cancer_type].shape[0]
            fig.add_annotation(
                x=cancer_type,
                y=count + 1,  
                text=str(count),
                showarrow=False,
                font=dict(size=10),
            )

        figs.append(fig)
    plot_htmls = [fig.to_html(full_html=False) for fig in figs]

    return render_template("histo.html", plot_htmls=plot_htmls)

if __name__ == "__main__":
    app.run(debug=True)