// Load environment variables from the google_key.env file
require('dotenv').config({ path: 'google_key.env' });

// Access the keys
const apiKey1 = process.env.GOOGLE_API_KEY_1;
const apiKey2 = process.env.GOOGLE_API_KEY_2;
const apiKey3 = process.env.GOOGLE_API_KEY_3;
const apiKey4 = process.env.GOOGLE_API_KEY_4;
const apiKey5 = process.env.GOOGLE_API_KEY_5;
const apiKey6 = process.env.GOOGLE_API_KEY_6;
const apiKey7 = process.env.GOOGLE_API_KEY_7;
const apiKey8 = process.env.GOOGLE_API_KEY_8;
const apiKey9 = process.env.GOOGLE_API_KEY_9;
const apiKey10 = process.env.GOOGLE_API_KEY_10;
const apiKey11 = process.env.GOOGLE_API_KEY_11;
const apiKey12 = process.env.GOOGLE_API_KEY_12;
const apiKey13 = process.env.GOOGLE_API_KEY_13;
const apiKey14 = process.env.GOOGLE_API_KEY_14;
const apiKey15 = process.env.GOOGLE_API_KEY_15;
const apiKey16 = process.env.GOOGLE_API_KEY_16;
const apiKey17 = process.env.GOOGLE_API_KEY_17;


// Function to choose an API key (this can be based on any logic you prefer)
function chooseApiKey() {
    // Example: Randomly choose one of the API keys
    const apiKeys = [apiKey1, apiKey2, apiKey3, apiKey4, apiKey5, apiKey6, apiKey7, apiKey8, apiKey9, apiKey10, apiKey11, apiKey12, apiKey13, apiKey14, apiKey15, apiKey16, apiKey17];
    const randomIndex = Math.floor(Math.random() * apiKeys.length);
    return apiKeys[randomIndex];
}

// Use the chosen API key
const selectedApiKey = chooseApiKey();
console.log('Selected API Key:', selectedApiKey);

// Example of using the selected API key with your existing configuration
generation_config = { 
    "max_output_tokens": 8192,
    "response_mime_type": "application/json",
    "api_key": selectedApiKey // Add the selected API key to the configuration
};