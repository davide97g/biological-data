
# ? Plot enriched terms in a word cloud.
from wordcloud import WordCloud
import matplotlib.pyplot as plt


def create_word_cloud(enriched_terms):

    wordcloud = WordCloud(background_color="white", width=1000, height=1000,
                          relative_scaling=0.5, normalize_plurals=False).generate_from_frequencies(enriched_terms)
    plt.figure(figsize=(8, 8), facecolor=None)
    plt.imshow(wordcloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.savefig("../../data/function/wordcloud.png")
    plt.show()
