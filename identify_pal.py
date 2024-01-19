import streamlit as st

def is_palindrome(s):
    return s == s[::-1]

def find_palindromes(input_string):
    palindromes = []

    # Iterate through all possible substrings
    for i in range(len(input_string)):
        for j in range(i + 40, len(input_string) + 40, 40):
            substring = input_string[i:j]

            # Check if the substring is a palindrome
            if is_palindrome(substring):
                palindromes.append(substring)

    return palindromes

# Streamlit application
# Page title
st.markdown(
    """
    # Identify palindromic sequences
    """
)

with st.form("ml-form"):
    file = st.file_uploader("FILE UPLOADER: Input the FASTA file of the DNA of your target", type='.fasta')
    st.markdown("**OR**")
    seq = st.text_area("Input your target DNA sequence")
    submitted = st.form_submit_button("Submit!")
